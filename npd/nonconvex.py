import logging
from fractions import Fraction
from itertools import chain
from typing import (Dict,
                    Tuple)

from gon.base import (Contour,
                      Point,
                      Segment)

from draw import draw
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from npd.utils import (fast_length,
                       inverse_shoelace)

LOGGER = logging.getLogger(__name__)

PLOT_STARTING_NUMBER = 20


def split(partition: Partition,
          *,
          requirement: Fraction) -> Tuple[Chunk, Chunk]:
    while not is_satisfactory_chunk_found(partition, requirement=requirement):
        join_pairs(partition)
        min_chunk_size = min(len(chunk.triangles)
                             for chunk in partition.chunk_map)
        LOGGER.debug(f"Min chunk size: {min_chunk_size}")
        if min_chunk_size >= 5:
            swap_leaves(partition)
    LOGGER.debug("Satisfactory chunk found!")
    join_small_chunks(partition, requirement=requirement)
    LOGGER.debug("Small chunks are joined")
    chunk, remainder = most_compact_pair(partition)
    if chunk.area > requirement:
        LOGGER.debug("chunk.area > requirement")
        chunk, remainder = readjust(chunk,
                                    remainder,
                                    requirement=requirement,
                                    partition=partition)
    elif chunk.area < requirement:
        LOGGER.debug("chunk.area < requirement")
        chunk, remainder = add_triangles(chunk,
                                         remainder,
                                         requirement=requirement,
                                         partition=partition)
    return chunk, remainder


def is_satisfactory_chunk_found(partition: Partition,
                                *,
                                requirement: Fraction) -> bool:
    """
    Checks if there is at least one chunk in the partition that
    has a satisfactory area for a given area requirement
    :param partition: input partition
    :param requirement: input area requirement
    :return: True if there is at least one satisfactory chunk
    """
    return any(is_satisfactory_area(chunk.area, requirement=requirement)
               for chunk in partition.chunk_map)


def is_satisfactory_area(area: Fraction,
                         *,
                         requirement) -> bool:
    """
    Checks if the area is close enough to the area requirement
    :param area: area of a polygon
    :param requirement: area requirement
    :return: True if the areas area close enough
    """
    return requirement - area <= 2 * area - requirement >= 0


def join_pairs(partition: Partition) -> None:
    """
    :param partition:
    :return:
    """
    smallest_chunk_size = min(len(chunk.triangles)
                              for chunk in partition.chunk_map)
    potential_sources = dict.fromkeys(
        chunk for chunk in partition.chunk_map
        if len(chunk.triangles) == smallest_chunk_size)
    potential_targets = dict.fromkeys(partition.chunk_map)
    while potential_sources:
        merge_least_compact(partition,
                            potential_sources,
                            potential_targets)


def merge_least_compact(partition: Partition,
                        potential_sources: Dict[Chunk, None],
                        potential_targets: Dict[Chunk, None]) -> None:
    """
    :param partition:
    :param potential_sources:
    :param potential_targets:
    :return:
    """
    source = min(potential_sources, key=Chunk.compactness)
    best_target = None
    max_compactness = 0
    neighbors = [neighbor for neighbor in partition.chunk_neighbors(source)
                 if neighbor in potential_targets]
    if not neighbors:
        potential_sources.pop(source)
        return
    for neighbor in neighbors:
        area = source.area + neighbor.area
        perimeter = partition.union_perimeter(source, neighbor)
        compactness_ = area ** 0.5 / perimeter
        if compactness_ > max_compactness:
            best_target = neighbor
            max_compactness = compactness_
    potential_sources.pop(source)
    potential_sources.pop(best_target, None)
    potential_targets.pop(source)
    potential_targets.pop(best_target)
    if len(source.triangles) > PLOT_STARTING_NUMBER:
        draw(*partition.chunk_map, source, best_target)
    partition.unite(source, best_target)


def swap_leaves(partition: Partition) -> None:
    """
    Locates leaf triangles in chunks and moves them to the neighbor
    chunk if this triangle is surrounded by two triangles from that
    chunk
    :param partition: input partition that will be updated in-place
    """
    for chunk in partition.chunk_map:
        triangle_to_move = None
        neighbor_chunk = None
        for triangle in chunk.triangles:
            triangle_neighbors = [
                neighbor
                for neighbor in partition.triangle_neighbors(triangle)]
            if len(triangle_neighbors) != 3:
                continue
            if sum(neighbor in chunk.triangles
                   for neighbor in triangle_neighbors) != 1:
                continue
            first_neighbor, second_neighbor = (
                neighbor for neighbor in triangle_neighbors
                if neighbor not in chunk.triangles)
            chunk_neighbors = partition.chunk_neighbors(chunk)
            first_neighbor_chunk = next(
                chunk_neighbor for chunk_neighbor in chunk_neighbors
                if first_neighbor in chunk_neighbor.triangles)
            second_neighbor_chunk = next(
                chunk_neighbor for chunk_neighbor in chunk_neighbors
                if second_neighbor in chunk_neighbor.triangles)
            if first_neighbor_chunk != second_neighbor_chunk:
                continue
            triangle_to_move = triangle
            neighbor_chunk = first_neighbor_chunk
            break
        if triangle_to_move is not None:
            partition.move(triangle_to_move,
                           source=chunk,
                           target=neighbor_chunk)


def join_small_chunks(partition: Partition,
                      *,
                      requirement: Fraction) -> None:
    """
    Finds chunks that are not close enough to the area requirement
    and joins them with their neighbors to form larger chunks.
    :param partition: input partition
    :param requirement: area requirement
    """
    potential_sources = dict.fromkeys(
        chunk for chunk in partition.chunk_map
        if not is_satisfactory_area(chunk.area, requirement=requirement))
    while potential_sources and len(partition.chunk_map) > 2:
        source = min(potential_sources, key=Chunk.compactness)
        best_target = None
        max_compactness = 0
        neighbors = list(partition.chunk_neighbors(source))
        for neighbor in neighbors:
            area = source.area + neighbor.area
            if area - requirement > requirement - source.area:
                continue
            perimeter = partition.union_perimeter(source, neighbor)
            compactness_ = area ** 0.5 / perimeter
            if compactness_ > max_compactness:
                best_target = neighbor
                max_compactness = compactness_
        potential_sources.pop(source)
        potential_sources.pop(best_target, None)
        if best_target is not None:
            union = partition.unite(source, best_target)
            if (not requirement - union.area <= 2 * union.area - requirement
                    >= 0):
                potential_sources.setdefault(union)


def most_compact_pair(partition: Partition) -> Tuple[Chunk, Chunk]:
    """
    For the given partition, finds the most compact pair of chunks
    :param partition: input partition
    :return: a pair of the most compact chunks
    """
    perimeter_edges = set(
        edge for triangle, neighbor_map in partition.triangular_map.items()
        for edge, neighbor in neighbor_map.items()
        if neighbor is None)
    total_perimeter = sum(map(fast_length, perimeter_edges))
    max_compactness = 0
    best_chunk = None
    for chunk in partition.chunk_map:
        if partition.is_articulation_chunk(chunk):
            continue
        chunk_compactness = chunk.area ** 0.5 / chunk.perimeter
        remainder_area = partition.total_area - chunk.area
        chunk_edges = {edge for triangle in chunk.triangles
                       for edge in partition.triangular_map[triangle]}
        common_edges = set(edge for edge in chunk_edges
                           if edge in perimeter_edges)
        chunk_only_edges = set(
            edge
            for triangle in chunk.triangles
            for edge, neighbor in partition.triangular_map[triangle].items()
            if neighbor is not None and neighbor not in chunk.triangles)
        chunk_only_perimeter = sum(map(fast_length, chunk_only_edges))
        common_perimeter = sum(map(fast_length, common_edges))
        remainder_perimeter = (total_perimeter - common_perimeter
                               + chunk_only_perimeter)
        remainder_compactness = remainder_area ** 0.5 / remainder_perimeter
        compactness_min_of_two = min(chunk_compactness, remainder_compactness)
        draw(*partition.chunk_map, red=chunk)
        if compactness_min_of_two > max_compactness:
            best_chunk = chunk
            max_compactness = compactness_min_of_two
    if best_chunk is None:
        raise RuntimeError("Couldn't find the best chunk")
    remainder = collect_remainder(best_chunk, partition)
    return best_chunk, remainder


def collect_remainder(chunk: Chunk,
                      partition: Partition) -> Chunk:
    """
    :param chunk:
    :param partition:
    :return:
    """
    while len(partition.chunk_map) > 2:
        neighbor = partition.chunk_neighbors(chunk)[0]
        neighbors_neighbor = next(
            node for node in partition.chunk_neighbors(neighbor)
            if node != chunk)
        partition.unite(neighbor, neighbors_neighbor)
    return partition.chunk_neighbors(chunk)[0]


def readjust(chunk: Chunk,
             neighbor: Chunk,
             *,
             requirement: Fraction,
             partition: Partition) -> Tuple[Chunk, Chunk]:
    """
    :param chunk:
    :param neighbor:
    :param requirement:
    :param partition:
    :return:
    """
    draw(*partition.chunk_map, red=chunk)
    ignored_neighbors = {}
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
            # border = dict.fromkeys(partition.border_triangles(chunk))
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
        elif chunk.area < requirement:
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
        if neighbor_triangles_count != 1 or chunk_triangles_count != 2:
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


def divide(partition: Partition,
           node: CachedPolygon,
           pivot_index: int,
           new_vertex: Point,
           chunk: Chunk,
           neighbor: Chunk) -> Tuple[CachedPolygon, CachedPolygon]:
    """
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
    left, right = divide(partition,
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
        return triangle.border.vertices.index(oriented_edge.end)
    neighbor = next(
        (neighbor for neighbor in partition.triangle_neighbors(triangle)
         if neighbor not in chunk.triangles), None)
    if neighbor is None:
        return 0  # any index will do
    unoriented_edge = partition.edge_between_triangles(triangle, neighbor)
    oriented_edge = next(edge for edge in triangle.edges
                         if edge == unoriented_edge)
    return triangle.border.vertices.index(oriented_edge.start)


def add_triangles(chunk: Chunk,
                  neighbor: Chunk,
                  *,
                  requirement: Fraction,
                  partition: Partition) -> Tuple[Chunk, Chunk]:
    """
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
