"""Splitting of non-convex polygons by merging triangles into chunks."""
import logging
from fractions import Fraction
from itertools import chain

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
        # Set-aside triangles can have been moved out of the chunk since.
        prioritized_chunk_triangles = chain(
            unignored_chunk_triangles,
            (triangle for triangle in ignored_neighbors
             if triangle in chunk.triangles))
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
           chunk=chunk)


def border_edge_triangle(chunk: Chunk,
                         neighbor: Chunk,
                         partition: Partition) -> CachedPolygon:
    """
    Finds a triangle of the chunk to split so that moves can continue.

    Used when no triangle can be moved between the chunks without
    disconnecting one of them. Returns the first triangle T of the chunk
    such that:
    1) T has exactly one neighbor triangle N in the neighbor chunk and
    the other two in the chunk itself (so no edge of T is on the
    polygon's border);
    2) at most one triangle of the neighbor chunk other than N shares a
    vertex with T and touches a chunk other than the neighbor (that is
    the given chunk: there are only two chunks when this is called).
    :param chunk: chunk to find the triangle in
    :param neighbor: neighbor chunk
    :param partition: partition the chunks belong to
    :return: the first triangle of the chunk meeting both conditions
    :raises ValueError: if no triangle of the chunk meets them
    """
    for triangle in chunk.triangles:
        triangle_neighbors = partition.triangle_neighbors(triangle)
        outside_neighbors = [triangle_neighbor
                             for triangle_neighbor in triangle_neighbors
                             if triangle_neighbor in neighbor.triangles]
        inside_neighbors_count = sum(
            triangle_neighbor in chunk.triangles
            for triangle_neighbor in triangle_neighbors)
        if (len(outside_neighbors) != 1
                or inside_neighbors_count != TRIANGLE_SIDES_COUNT - 1):
            continue
        (outside_neighbor,) = outside_neighbors
        vertices = triangle.border.vertices
        touching_triangles = [
            neighbor_triangle for neighbor_triangle in neighbor.triangles
            if neighbor_triangle is not outside_neighbor
            and any(vertex in vertices
                    for vertex in neighbor_triangle.border.vertices)
            and partition.is_touching_other_chunk(neighbor_triangle,
                                                  this_chunk=neighbor)]
        draw(*partition.chunk_map, outside_neighbor,
             blue=chunk, red=triangle)
        if len(touching_triangles) <= 1:
            return triangle
    raise ValueError("Couldn't find a triangle to split")


def edge_neighbor(triangle: CachedPolygon,
                  chunk: Chunk,
                  partition: Partition) -> CachedPolygon:
    """
    Finds the neighbor of the triangle to split along with it.

    Prefers the first neighbor C of the triangle within the chunk such
    that no neighbor of C other than the triangle touches another chunk,
    i.e. C reaches the chunk's border only through the triangle; without
    one, takes the first neighbor within the chunk. Either works for
    ``divide_border_edge_triangle``: cutting the triangle at the midpoint
    of its edge with any neighbor in the chunk leaves a half that
    touches the other chunk and can be moved without disconnecting
    anything (the other half and the neighbor's two halves keep the rest
    connected). Expects a triangle touching another chunk, as
    ``border_edge_triangle`` returns.
    :param triangle: triangle of the chunk touching another chunk
    :param chunk: chunk the triangle belongs to
    :param partition: partition the chunk belongs to
    :return: the first neighbor of the triangle within the chunk meeting
    the condition, or the first one within the chunk if none does
    :raises ValueError: if the triangle has no neighbor within the chunk
    """
    chunk_neighbors = [neighbor
                       for neighbor in partition.triangle_neighbors(triangle)
                       if neighbor in chunk.triangles]
    if not chunk_neighbors:
        raise ValueError("Couldn't find edge neighbor.")
    for neighbor in chunk_neighbors:
        neighbor_neighbors = partition.triangle_neighbors(neighbor)
        # `triangle` is among them and touches another chunk, so a count
        # of 1 means that none of the others does.
        count = sum(partition.is_touching_other_chunk(neighbor_neighbor,
                                                      this_chunk=chunk)
                    for neighbor_neighbor in neighbor_neighbors)

        draw(*partition.chunk_map, blue=chunk, red=neighbor)
        if count == 1:
            return neighbor
    return chunk_neighbors[0]


def divide(partition: Partition,
           *,
           node: CachedPolygon,
           pivot_index: int,
           new_vertex: Point[Fraction],
           chunk: Chunk) -> tuple[CachedPolygon, CachedPolygon]:
    """
    Cuts the node by a segment from its pivot vertex to the new vertex.

    The new vertex must lie strictly inside the node's edge opposite the
    pivot. The triangle across that edge, if any, is cut at the new
    vertex too, so that no vertex lies inside an edge. Updates the
    triangular map, the chunks' triangles and the edges chunks share.
    Areas don't change; neither do the cached perimeters, though cutting
    an edge on the polygon's border can change their rounding.
    :param partition: partition to update in place
    :param node: triangle of the chunk to cut
    :param pivot_index: index of the pivot among the node's vertices
    :param new_vertex: point strictly inside the edge opposite the pivot
    :param chunk: chunk the node belongs to
    :return: the part with the vertex before the pivot and the part with
    the vertex after it
    :raises ValueError: if the node doesn't belong to the chunk
    """
    if node not in chunk.triangles:
        raise ValueError("The node must belong to the chunk")
    triangular_map = partition.triangular_map
    vertices = node.border.vertices
    pivot_vertex = vertices[pivot_index]
    left_vertex = vertices[pivot_index - 1]
    right_vertex = vertices[pivot_index - 2]
    left_segment = Segment(pivot_vertex, left_vertex)
    right_segment = Segment(right_vertex, pivot_vertex)
    third_segment = Segment(left_vertex, right_vertex)
    left_third_segment = Segment(left_vertex, new_vertex)
    right_third_segment = Segment(right_vertex, new_vertex)
    cut = Segment(pivot_vertex, new_vertex)
    left_part = CachedPolygon(Contour([pivot_vertex, left_vertex, new_vertex]))
    right_part = CachedPolygon(Contour([right_vertex, pivot_vertex,
                                        new_vertex]))
    node_neighbors = triangular_map.pop(node)
    triangular_map[left_part] = {cut: right_part}
    triangular_map[right_part] = {cut: left_part}
    _link(partition, left_part, node_neighbors[left_segment], left_segment)
    _link(partition, right_part, node_neighbors[right_segment],
          right_segment)
    third_neighbor = node_neighbors[third_segment]
    if third_neighbor is None:
        triangular_map[left_part][left_third_segment] = None
        triangular_map[right_part][right_third_segment] = None
    else:
        only_vertex = next(vertex
                           for vertex in third_neighbor.border.vertices
                           if vertex not in {left_vertex, right_vertex})
        third_cut = Segment(new_vertex, only_vertex)
        third_left_segment = Segment(left_vertex, only_vertex)
        third_right_segment = Segment(right_vertex, only_vertex)
        third_left_part = CachedPolygon(
            Contour([new_vertex, left_vertex, only_vertex]))
        third_right_part = CachedPolygon(
            Contour([right_vertex, new_vertex, only_vertex]))
        third_neighbor_neighbors = triangular_map.pop(third_neighbor)
        triangular_map[third_left_part] = {third_cut: third_right_part}
        triangular_map[third_right_part] = {third_cut: third_left_part}
        _link(partition, left_part, third_left_part, left_third_segment)
        _link(partition, right_part, third_right_part, right_third_segment)
        _link(partition, third_left_part,
              third_neighbor_neighbors[third_left_segment],
              third_left_segment)
        _link(partition, third_right_part,
              third_neighbor_neighbors[third_right_segment],
              third_right_segment)
        owner = next(chunk_ for chunk_ in partition.chunk_map
                     if third_neighbor in chunk_.triangles)
        owner.triangles.pop(third_neighbor)
        owner.triangles.setdefault(third_left_part)
        owner.triangles.setdefault(third_right_part)
        if owner is not chunk:
            partition.split_shared_edge(chunk,
                                        owner,
                                        third_segment,
                                        [left_third_segment,
                                         right_third_segment])
    chunk.triangles.pop(node)
    chunk.triangles.setdefault(left_part)
    chunk.triangles.setdefault(right_part)
    return left_part, right_part


def _link(partition: Partition,
          triangle: CachedPolygon,
          other: CachedPolygon | None,
          edge: Segment[Fraction]) -> None:
    """Records what lies across the edge, on both sides if needed."""
    partition.triangular_map[triangle][edge] = other
    if other is not None:
        partition.triangular_map[other][edge] = triangle


def final_node_refinement(node_to_split: CachedPolygon,
                          *,
                          requirement: Fraction,
                          chunk: Chunk,
                          neighbor: Chunk,
                          partition: Partition) -> None:
    """
    Cuts the node so that the chunk's area meets the requirement exactly.

    Cuts the node from the pivot ``find_pivot`` picks, so that the part
    staying in the chunk has the area the chunk lacks without the node,
    and moves the other part to the neighbor. Callers have just moved
    the node into the chunk, past the requirement: ``readjust`` moved it
    out, fell below the requirement and moved it back; ``add_triangles``
    moved it in and rose above it. The chunk stays connected if the node
    either has a neighbor in it and isn't an articulation triangle, or
    joined it from outside; the moved part touches the neighbor if the
    node does.
    :param node_to_split: triangle of the chunk to cut
    :param requirement: required area of the chunk
    :param chunk: chunk the triangle belongs to, updated in place
    :param neighbor: chunk touching the triangle, receiving the cut-off
    part
    :param partition: partition the chunks belong to
    :raises ValueError: if the node isn't in the chunk, or the
    requirement isn't strictly between the chunk's area without the node
    and with it
    """
    if node_to_split not in chunk.triangles:
        raise ValueError("The node must belong to the chunk")
    if not chunk.area - node_to_split.area < requirement < chunk.area:
        raise ValueError("The requirement must be strictly between the "
                         "chunk's area without the node and with it")
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
                      chunk=chunk)
    partition.move(right,
                   source=chunk,
                   target=neighbor)


def find_pivot(triangle: CachedPolygon,
               *,
               partition: Partition,
               chunk: Chunk) -> int:
    """
    Finds the index of the vertex to cut the triangle from.

    ``final_node_refinement`` cuts the triangle from the pivot P to the
    edge opposite it, keeps the part with the vertex before P in the
    chunk and moves the part with the vertex after P to the other chunk.
    So that the kept part stays attached to the chunk, the edge from the
    vertex before P to P is the one shared with the triangle's first
    neighbor in the chunk. If there is none, the moved part stays
    attached to the other chunk instead: the edge from P to the vertex
    after it is the one shared with the triangle's first neighbor
    outside the chunk (the other chunk, as there are only two chunks
    when this is called).
    :param triangle: triangle of the chunk to cut
    :param partition: partition the chunk belongs to
    :param chunk: chunk the triangle belongs to
    :return: index of the pivot among the triangle's vertices, 0 if the
    triangle has no neighbors
    """
    neighbors = partition.triangular_map[triangle]
    vertices = triangle.border.vertices
    chunk_edge = next((edge for edge, neighbor in neighbors.items()
                       if neighbor is not None
                       and neighbor in chunk.triangles), None)
    if chunk_edge is not None:
        return next(index for index in range(len(vertices))
                    if Segment(vertices[index - 1], vertices[index])
                    == chunk_edge)
    outside_edge = next((edge for edge, neighbor in neighbors.items()
                         if neighbor is not None), None)
    if outside_edge is None:
        return 0  # any index will do
    return next(index for index in range(len(vertices))
                if Segment(vertices[index], vertices[index - 2])
                == outside_edge)


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
