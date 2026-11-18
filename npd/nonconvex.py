"""Splitting of non-convex polygons by merging triangles into chunks."""
import logging
from fractions import Fraction
from itertools import chain
from typing import cast

from gon.base import (Contour,
                      Point,
                      Segment)

from npd.draw import draw
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from npd.utils import (fast_length,
                       inverse_shoelace)

LOGGER = logging.getLogger(__name__)

PLOT_STARTING_NUMBER = 20
# leaves are swapped only once all chunks have at least this many triangles
LEAF_SWAPPING_MIN_CHUNK_SIZE = 5
# the partition is reduced to a chunk and the remainder
PARTS_COUNT = 2
TRIANGLE_SIDES_COUNT = 3


def split(partition: Partition,
          *,
          requirement: Fraction) -> tuple[Chunk, Chunk]:
    """
    Splits the partition into a chunk of the required area and the rest.

    :param partition: partition of the polygon, updated in-place
    :param requirement: area of the chunk to split off
    :return: the chunk with the required area and the remainder chunk
    :raises ValueError: if the requirement isn't positive and less than
    the partition's total area
    """
    if not 0 < requirement < partition.total_area:
        raise ValueError("Requirement must be positive and less than "
                         "the partition's total area.")
    # Joining stops at two chunks, even if neither is satisfactory (only
    # possible for requirements above 3/4 of the total area): the steps
    # below need a chunk and a remainder.
    while (len(partition.chunk_map) > PARTS_COUNT
           and not is_satisfactory_chunk_found(partition,
                                               requirement=requirement)):
        join_pairs(partition)
        min_chunk_size = min(len(chunk.triangles)
                             for chunk in partition.chunk_map)
        if min_chunk_size >= LEAF_SWAPPING_MIN_CHUNK_SIZE:
            swap_leaves(partition)
    join_small_chunks(partition, requirement=requirement)
    chunk, remainder = most_compact_pair(partition)
    if chunk.area > requirement:
        chunk, remainder = readjust(chunk,
                                    remainder,
                                    requirement=requirement,
                                    partition=partition)
    elif chunk.area < requirement:
        chunk, remainder = add_triangles(chunk,
                                         remainder,
                                         requirement=requirement,
                                         partition=partition)
    return chunk, remainder


def is_satisfactory_chunk_found(partition: Partition,
                                *,
                                requirement: Fraction) -> bool:
    """
    Checks if any chunk is large enough to stop joining chunks.

    See ``is_satisfactory_area`` for when a chunk is large enough.

    :param partition: input partition
    :param requirement: input area requirement
    :return: True if there is at least one satisfactory chunk
    """
    return any(is_satisfactory_area(chunk.area, requirement=requirement)
               for chunk in partition.chunk_map)


def is_satisfactory_area(area: Fraction,
                         *,
                         requirement: Fraction) -> bool:
    """
    Checks if chunks should stop being joined to approach the requirement.

    Joining pairs of chunks roughly doubles their area. It's worth
    stopping once the next doubling would overshoot the requirement and
    leave a gap at least as large as the current one::

        |requirement - area| <= |2 * area - requirement|
        and 2 * area >= requirement

    If ``area >= requirement``, this always holds, since doubling only
    moves further away. If ``requirement / 2 <= area < requirement``, it
    becomes ``requirement - area <= 2 * area - requirement``, that is
    ``area >= 2/3 * requirement``. If ``area < requirement / 2``, the next
    doubling still undershoots and gets closer, so it never holds. All
    cases together reduce to ``3 * area >= 2 * requirement``.

    The doubling is only approximate, as joined chunks can differ in
    size, so the 2/3 threshold is a heuristic, not a guarantee.

    :param area: area of a chunk
    :param requirement: area requirement
    :return: True if the chunk is close enough to stop joining
    """
    return 3 * area >= 2 * requirement


def join_pairs(partition: Partition) -> None:
    """
    Merges the smallest chunks with their neighbors into compact pairs.

    Goes through the chunks with the fewest triangles from the least
    compact one, uniting each with the neighbor that gives the most
    compact union. A chunk takes part in at most one union, and unions
    aren't merged further.
    :param partition: partition to update in place
    """
    smallest_chunk_size = min(len(chunk.triangles)
                              for chunk in partition.chunk_map)
    # Sorted once: merged chunks are skipped and unions never become
    # sources, so a source's compactness can't change while waiting.
    sources = sorted((chunk for chunk in partition.chunk_map
                      if len(chunk.triangles) == smallest_chunk_size),
                     key=Chunk.compactness)
    unmerged = dict.fromkeys(partition.chunk_map)
    for source in sources:
        if source not in unmerged:
            continue
        neighbors = [neighbor
                     for neighbor in partition.chunk_neighbors(source)
                     if neighbor in unmerged]
        if not neighbors:
            continue
        best_target = max(
            neighbors,
            key=lambda neighbor: partition.union_compactness(source,
                                                             neighbor))
        del unmerged[source]
        del unmerged[best_target]
        if len(source.triangles) > PLOT_STARTING_NUMBER:
            draw(*partition.chunk_map, source, best_target)
        partition.unite(source, best_target)


def swap_leaves(partition: Partition) -> None:
    """
    Moves leaf triangles to the neighbor chunks surrounding them.

    A leaf is a triangle with no edge on the polygon's border, one
    neighbor in its own chunk and the other two in a single other chunk.
    Moving it there keeps its chunk non-empty and connected, and
    shortens the perimeters of both chunks: the edge to its own chunk is
    shorter than the two edges to the other one. Triangles on the
    polygon's border are never moved: the target's perimeter would grow.
    At most one leaf, the first found, is moved out of each chunk per
    call. Chunks are visited in the order of the chunk map, each seeing
    the moves made before it.
    :param partition: input partition that will be updated in-place
    """
    for chunk in partition.chunk_map:
        for triangle in chunk.triangles:
            triangle_neighbors = partition.triangle_neighbors(triangle)
            if len(triangle_neighbors) != TRIANGLE_SIDES_COUNT:
                continue
            outside_neighbors = [neighbor
                                 for neighbor in triangle_neighbors
                                 if neighbor not in chunk.triangles]
            if len(outside_neighbors) != TRIANGLE_SIDES_COUNT - 1:
                continue
            first_neighbor, second_neighbor = outside_neighbors
            chunk_neighbors = partition.chunk_neighbors(chunk)
            first_neighbor_chunk = next(
                chunk_neighbor for chunk_neighbor in chunk_neighbors
                if first_neighbor in chunk_neighbor.triangles)
            second_neighbor_chunk = next(
                chunk_neighbor for chunk_neighbor in chunk_neighbors
                if second_neighbor in chunk_neighbor.triangles)
            if first_neighbor_chunk is not second_neighbor_chunk:
                continue
            # safe to mutate `chunk.triangles`: we stop iterating right after
            partition.move(triangle,
                           source=chunk,
                           target=first_neighbor_chunk)
            break


def join_small_chunks(partition: Partition,
                      *,
                      requirement: Fraction) -> None:
    """
    Joins chunks with unsatisfactory areas with their neighbors.

    Repeatedly takes the least compact chunk whose area isn't
    satisfactory (see ``is_satisfactory_area``) and unites it with the
    neighbor giving the most compact union. Neighbors whose union would
    exceed the requirement by more than the chunk falls short of it are
    skipped, and a chunk with no other neighbor is left as is. A union
    that is still unsatisfactory is joined further. Stops when no
    unsatisfactory chunk can be joined or only two chunks are left.
    :param partition: partition to update in place
    :param requirement: area requirement
    """
    sources = dict.fromkeys(
        chunk for chunk in partition.chunk_map
        if not is_satisfactory_area(chunk.area, requirement=requirement))
    while sources and len(partition.chunk_map) > PARTS_COUNT:
        source = min(sources, key=Chunk.compactness)
        del sources[source]
        shortfall = requirement - source.area
        targets = [neighbor
                   for neighbor in partition.chunk_neighbors(source)
                   if source.area + neighbor.area - requirement <= shortfall]
        if not targets:
            continue
        best_target = max(
            targets,
            key=lambda neighbor: partition.union_compactness(source,
                                                             neighbor))
        sources.pop(best_target, None)
        union = partition.unite(source, best_target)
        if not is_satisfactory_area(union.area, requirement=requirement):
            sources.setdefault(union)


def most_compact_pair(partition: Partition) -> tuple[Chunk, Chunk]:
    """
    Picks a chunk and unites all the other chunks into its remainder.

    Among the chunks whose removal leaves the other chunks connected,
    picks the one for which the less compact of the chunk and its
    remainder is the most compact (the first one on ties).
    :param partition: partition with at least two chunks, updated
    in-place: all the chunks except the picked one are united
    :return: the picked chunk and the remainder
    """
    border_perimeter = sum(fast_length(edge)
                           for neighbors in partition.triangular_map.values()
                           for edge, neighbor in neighbors.items()
                           if neighbor is None)

    def score(chunk: Chunk) -> float:
        # Each edge appears once among the chunk's triangles, so there's
        # no need to collect them in sets (hashing them is slow).
        chunk_border_perimeter = shared_perimeter = 0.
        for triangle in chunk.triangles:
            for edge, neighbor in partition.triangular_map[triangle].items():
                if neighbor is None:
                    chunk_border_perimeter += fast_length(edge)
                elif neighbor not in chunk.triangles:
                    shared_perimeter += fast_length(edge)
        # The remainder's border is the polygon's border without the
        # chunk's part of it, plus the edges it shares with the chunk.
        remainder_perimeter = (border_perimeter - chunk_border_perimeter
                               + shared_perimeter)
        remainder_compactness: float = (
            (partition.total_area - chunk.area) ** 0.5
            / remainder_perimeter)
        return min(chunk.compactness(), remainder_compactness)

    best_chunk = max((chunk for chunk in partition.chunk_map
                      if not partition.is_articulation_chunk(chunk)),
                     key=score)
    draw(*partition.chunk_map, red=best_chunk)
    return best_chunk, collect_remainder(best_chunk, partition)


def collect_remainder(chunk: Chunk,
                      partition: Partition) -> Chunk:
    """
    Unites all the chunks except the given one into a single chunk.

    Repeatedly unites the chunk's first neighbor with one of that
    neighbor's own neighbors until only the chunk and one other chunk
    are left.
    :param chunk: chunk to keep out of the union
    :param partition: partition with at least two chunks, updated
    in-place
    :return: the union of all the other chunks
    :raises ValueError: if the partition has fewer than two chunks or
    the other chunks aren't connected without the given one
    """
    if partition.is_articulation_chunk(chunk):
        raise ValueError("The other chunks must be connected without "
                         "the given one")
    while len(partition.chunk_map) > PARTS_COUNT:
        neighbor = partition.chunk_neighbors(chunk)[0]
        neighbors_neighbor = next(
            node for node in partition.chunk_neighbors(neighbor)
            if node is not chunk)
        partition.unite(neighbor, neighbors_neighbor)
    return partition.chunk_neighbors(chunk)[0]


def readjust(chunk: Chunk,
             neighbor: Chunk,
             *,
             requirement: Fraction,
             partition: Partition) -> tuple[Chunk, Chunk]:
    """
    Moves triangles from the chunk to its neighbor to meet the requirement.

    :param chunk:
    :param neighbor:
    :param requirement:
    :param partition:
    :return:
    """
    draw(*partition.chunk_map, red=chunk)
    ignored_neighbors: dict[CachedPolygon, None] = {}
    while True:
        unignored_chunk_triangles = (triangle for triangle in chunk.triangles
                                     if triangle not in ignored_neighbors)
        prioritized_chunk_triangles = chain(unignored_chunk_triangles,
                                            ignored_neighbors)
        border_triangle = next(
            (triangle for triangle in prioritized_chunk_triangles
             if partition.is_touching_other_chunk(triangle, this_chunk=chunk)
             and not partition.is_articulation_triangle(triangle,
                                                        chunk.triangles)
             ), None)
        if border_triangle is None and ignored_neighbors:
            ignored_neighbors = {}
            continue
        if border_triangle is None:
            divide_border_edge_triangle(chunk, neighbor, partition)
            border_triangle = next(
                triangle for triangle in partition.border_triangles(chunk)
                if not partition.is_articulation_triangle(triangle,
                                                          chunk.triangles))
            draw(*partition.chunk_map, blue=chunk)

        partition.move(border_triangle, source=chunk, target=neighbor)
        for neighbor_triangle in partition.triangle_neighbors(border_triangle):
            if neighbor_triangle in chunk.triangles:
                ignored_neighbors.setdefault(neighbor_triangle)
        draw(*partition.chunk_map, blue=chunk)
        if chunk.area == requirement:
            return chunk, neighbor
        if chunk.area < requirement:
            partition.move(border_triangle,
                           source=neighbor,
                           target=chunk)
            final_node_refinement(border_triangle,
                                  requirement=requirement,
                                  chunk=chunk,
                                  neighbor=neighbor,
                                  partition=partition)
            return chunk, neighbor


def divide_border_edge_triangle(chunk: Chunk,
                                neighbor: Chunk,
                                partition: Partition) -> None:
    """
    Splits a border-edge triangle of the chunk into two triangles.

    :param chunk: chunk where a node on the border will be found
    that will be split
    :param neighbor: neighbor chunk
    :param partition: partition
    """
    leaf = border_edge_triangle(chunk, neighbor, partition)
    neighbor_to_split = edge_neighbor(leaf, chunk, partition)
    edge = partition.edge_between_triangles(leaf, neighbor_to_split)
    pivot_index = next(index
                       for index, vertex in enumerate(leaf.border.vertices)
                       if vertex not in edge)
    divide(partition,
           node=leaf,
           pivot_index=pivot_index,
           new_vertex=edge.centroid,
           chunk=chunk,
           neighbor=neighbor)


def border_edge_triangle(chunk: Chunk,
                         neighbor: Chunk,
                         partition: Partition) -> CachedPolygon:
    """
    Finds a border-edge triangle in the chunk.

    We call "border-edge" triangles those triangles that can be split
    so that the algorithm would be able to find a way to continue
    reassigning triangles from one chunk to another.
    These border-edge triangles have the following properties:
    1) they have exactly 1 neighbor triangle N that belongs to neighbor
    chunk and 2 neighbor triangles that belong to the current chunk
    2) there is not more than 1 triangle sharing any vertex with N
    that is adjacent to any triangle from the current chunk
    """
    for triangle in chunk.triangles:
        triangle_neighbors = partition.triangle_neighbors(triangle)
        neighbor_triangles_count = sum(
            triangle_neighbor in neighbor.triangles
            for triangle_neighbor in triangle_neighbors)
        chunk_triangles_count = sum(
            triangle_neighbor in chunk.triangles
            for triangle_neighbor in triangle_neighbors)
        if (neighbor_triangles_count != 1
                or chunk_triangles_count != TRIANGLE_SIDES_COUNT - 1):
            continue
        triangle_neighbor_triangle = next(
            triangle_neighbor
            for triangle_neighbor in triangle_neighbors
            if triangle_neighbor in neighbor.triangles)
        vertices = triangle.border.vertices
        triangles_sharing_vertices = {
            neighbor_triangle for neighbor_triangle in neighbor.triangles
            if neighbor_triangle != triangle_neighbor_triangle
            and any(vertex in vertices
                    for vertex in neighbor_triangle.border.vertices)}
        adjacent_triangles = {
            neighbor_triangle
            for neighbor_triangle in triangles_sharing_vertices
            if partition.is_touching_other_chunk(neighbor_triangle,
                                                 this_chunk=neighbor)}
        draw(*partition.chunk_map, triangle_neighbor_triangle,
             blue=chunk, red=triangle)
        if len(adjacent_triangles) <= 1:
            return triangle
    raise ValueError("Couldn't find a triangle to split")


def edge_neighbor(triangle: CachedPolygon,
                  chunk: Chunk,
                  partition: Partition) -> CachedPolygon:
    """
    Finds the neighbor of the triangle that should be split.

    :param triangle:
    :param chunk:
    :param partition:
    :return:
    """
    for neighbor in partition.triangle_neighbors(triangle):
        if neighbor not in chunk.triangles:
            continue
        neighbor_neighbors = partition.triangle_neighbors(neighbor)
        count = sum(partition.is_touching_other_chunk(neighbor_neighbor,
                                                      this_chunk=chunk)
                    for neighbor_neighbor in neighbor_neighbors)

        draw(*partition.chunk_map, blue=chunk, red=neighbor)
        if count == 1:
            return neighbor
    raise ValueError("Couldn't find edge neighbor.")


def divide(partition: Partition,  # noqa: PLR0917, PLR0913
           node: CachedPolygon,
           pivot_index: int,
           new_vertex: Point,
           chunk: Chunk,
           neighbor: Chunk) -> tuple[CachedPolygon, CachedPolygon]:
    """
    Splits the node by a segment from its pivot vertex to the new vertex.

    :param partition:
    :param node:
    :param pivot_index:
    :param new_vertex:
    :param chunk:
    :param neighbor:
    :return:
    """
    vertices = node.border.vertices
    pivot_vertex = vertices[pivot_index]
    left_vertex = vertices[pivot_index - 1]
    right_vertex = vertices[pivot_index - 2]
    left_segment = Segment(pivot_vertex, left_vertex)
    right_segment = Segment(right_vertex, pivot_vertex)
    third_segment = Segment(left_vertex, right_vertex)
    left_part = CachedPolygon(Contour([pivot_vertex, left_vertex, new_vertex]))
    right_part = CachedPolygon(Contour([right_vertex, pivot_vertex,
                                        new_vertex]))
    left_neighbor = partition.triangular_map[node][left_segment]
    right_neighbor = partition.triangular_map[node][right_segment]
    third_neighbor = partition.triangular_map[node][third_segment]
    partition.triangular_map[left_part] = {Segment(pivot_vertex, new_vertex):
                                           right_part}
    partition.triangular_map[right_part] = {Segment(pivot_vertex, new_vertex):
                                            left_part}
    partition.triangular_map[left_part][left_segment] = left_neighbor
    if left_neighbor is not None:
        partition.triangular_map[left_neighbor][left_segment] = left_part
    partition.triangular_map[right_part][right_segment] = right_neighbor
    if right_neighbor is not None:
        partition.triangular_map[right_neighbor][right_segment] = right_part
    if third_neighbor is not None:
        third_neighbor_only_vertex = next(
            vertex for vertex in third_neighbor.border.vertices
            if vertex not in {left_vertex, right_vertex})
        third_neighbor_left_part = CachedPolygon(
            Contour([new_vertex, left_vertex, third_neighbor_only_vertex]))
        third_neighbor_right_part = CachedPolygon(
            Contour([right_vertex, new_vertex, third_neighbor_only_vertex]))
        partition.triangular_map[third_neighbor_left_part] = {
            Segment(new_vertex, third_neighbor_only_vertex):
                third_neighbor_right_part}
        partition.triangular_map[third_neighbor_right_part] = {
            Segment(new_vertex, third_neighbor_only_vertex):
                third_neighbor_left_part}
        partition.triangular_map[left_part][
            Segment(left_vertex, new_vertex)] = third_neighbor_left_part
        partition.triangular_map[third_neighbor_left_part][
            Segment(left_vertex, new_vertex)] = left_part
        partition.triangular_map[right_part][
            Segment(right_vertex, new_vertex)] = third_neighbor_right_part
        partition.triangular_map[third_neighbor_right_part][
            Segment(right_vertex, new_vertex)] = right_part
        third_neighbor_left_segment = Segment(left_vertex,
                                              third_neighbor_only_vertex)
        third_neighbor_right_segment = Segment(right_vertex,
                                               third_neighbor_only_vertex)
        third_neighbor_left_neighbor = next(
            (neighbor for neighbor in partition.triangle_neighbors(
                third_neighbor)
             if partition.edge_between_triangles(third_neighbor, neighbor)
             == third_neighbor_left_segment), None)
        third_neighbor_right_neighbor = next(
            (neighbor for neighbor in partition.triangle_neighbors(
                third_neighbor)
             if partition.edge_between_triangles(third_neighbor, neighbor)
             == third_neighbor_right_segment), None)
        partition.triangular_map[third_neighbor_left_part][
            third_neighbor_left_segment] = third_neighbor_left_neighbor
        if third_neighbor_left_neighbor is not None:
            partition.triangular_map[third_neighbor_left_neighbor][
                third_neighbor_left_segment] = third_neighbor_left_part
        partition.triangular_map[third_neighbor_right_part][
            third_neighbor_right_segment] = third_neighbor_right_neighbor
        if third_neighbor_right_neighbor is not None:
            partition.triangular_map[third_neighbor_right_neighbor][
                third_neighbor_right_segment] = third_neighbor_right_part
        partition.triangular_map.pop(third_neighbor)
        chunk_with_third_neighbor = next(chunk_ for chunk_ in (chunk, neighbor)
                                         if third_neighbor in chunk_.triangles)
        chunk_with_third_neighbor.triangles.pop(third_neighbor)
        chunk_with_third_neighbor.triangles.setdefault(
            third_neighbor_left_part)
        chunk_with_third_neighbor.triangles.setdefault(
            third_neighbor_right_part)

    partition.triangular_map.pop(node)
    chunk.triangles.pop(node)
    chunk.triangles.setdefault(left_part)
    chunk.triangles.setdefault(right_part)
    return left_part, right_part


def final_node_refinement(node_to_split: CachedPolygon,
                          *,
                          requirement: Fraction,
                          chunk: Chunk,
                          neighbor: Chunk,
                          partition: Partition) -> None:
    """
    Splits the node so that the chunk area meets the requirement exactly.

    :param node_to_split: triangle to be split
    :param requirement: area requirement for this triangle
    :param chunk: chunk the triangle belongs to
    :param neighbor: neighbor chunk touching the triangle
    :param partition: partition
    """
    pivot_index = find_pivot(node_to_split,
                             partition=partition,
                             chunk=chunk)
    vertices = node_to_split.border.vertices
    draw(*partition.chunk_map, red=[node_to_split, vertices[pivot_index]],
         blue=chunk)
    new_vertex = inverse_shoelace(
        requirement + node_to_split.area - chunk.area,
        endpoint=vertices[pivot_index],
        base_start=vertices[pivot_index - 1],
        base_end=vertices[pivot_index - 2])
    draw(*partition.chunk_map,
         red=[node_to_split, vertices[pivot_index]], blue=[chunk, new_vertex])
    _, right = divide(partition,
                      node=node_to_split,
                      pivot_index=pivot_index,
                      new_vertex=new_vertex,
                      chunk=chunk,
                      neighbor=neighbor)
    partition.move(right,
                   source=chunk,
                   target=neighbor)


def find_pivot(triangle: CachedPolygon,
               *,
               partition: Partition,
               chunk: Chunk) -> int:
    """
    Finds the index of the vertex to split the triangle from.

    :param triangle:
    :param partition:
    :param chunk:
    :return:
    """
    predecessor = next(
        (neighbor for neighbor in partition.triangle_neighbors(triangle)
         if neighbor in chunk.triangles), None)
    if predecessor is not None:
        unoriented_edge = partition.edge_between_triangles(triangle,
                                                           predecessor)
        oriented_edge = next(edge for edge in triangle.edges
                             if edge == unoriented_edge)
        return cast("int",
                    triangle.border.vertices.index(oriented_edge.end))
    neighbor = next(
        (neighbor for neighbor in partition.triangle_neighbors(triangle)
         if neighbor not in chunk.triangles), None)
    if neighbor is None:
        return 0  # any index will do
    unoriented_edge = partition.edge_between_triangles(triangle, neighbor)
    oriented_edge = next(edge for edge in triangle.edges
                         if edge == unoriented_edge)
    return cast("int", triangle.border.vertices.index(oriented_edge.start))


def add_triangles(chunk: Chunk,
                  neighbor: Chunk,
                  *,
                  requirement: Fraction,
                  partition: Partition) -> tuple[Chunk, Chunk]:
    """
    Moves triangles from the neighbor to the chunk to meet the requirement.

    :param chunk:
    :param neighbor:
    :param requirement:
    :param partition:
    :return:
    """
    draw(*partition.chunk_map, red=chunk)
    while True:
        triangle = next(partition.movable_triangles(neighbor), None)
        if triangle is None:
            divide_border_edge_triangle(neighbor, chunk, partition)
            continue
        partition.move(triangle, source=neighbor, target=chunk)
        if chunk.area < requirement:
            continue
        if chunk.area == requirement:
            return chunk, neighbor
        final_node_refinement(triangle,
                              requirement=requirement,
                              chunk=chunk,
                              neighbor=neighbor,
                              partition=partition)
        return chunk, neighbor
