# """Same but with attempt to create my own class"""
# import math
# import random
# from collections import defaultdict
# from fractions import Fraction
# from functools import reduce
# from itertools import product
# from operator import (attrgetter,
#                       or_)
# from typing import (Dict,
#                     FrozenSet,
#                     Iterable,
#                     Iterator,
#                     List,
#                     Optional,
#                     Sequence,
#                     Set,
#                     Tuple,
#                     TypeVar)
#
# import matplotlib.pyplot as plt
# import networkx as nx
# from gon.base import (Compound,
#                       Contour,
#                       EMPTY,
#                       Linear,
#                       Orientation,
#                       Point,
#                       Polygon,
#                       Segment,
#                       Triangulation)
# from gon.hints import Scalar
# from ground.base import get_context
#
# from draw import (draw,
#                   show)
#
# T = TypeVar('T')
#
# context = get_context()
#
# DRAW = True
# # DRAW = False
#
# PLOT_STARTING_NUMBER = 5
#
#
# class CachedPolygon(Polygon):
#     slots = '_hash',
#
#     def __init__(self,
#                  border: Contour[Scalar],
#                  holes: Optional[Sequence[Contour[Scalar]]] = None) -> None:
#         super().__init__(border, holes)
#         self._hash = super().__hash__()
#         self._perimeter = (fast_length(self.border)
#                            + sum(map(fast_length, holes)))
#
#     def __hash__(self) -> int:
#         return self._hash
#
#     @property
#     def perimeter(self):
#         return self._perimeter
#
#
# class Chunk:
#     def __init__(self, *,
#                  triangles: Set[CachedPolygon],
#                  area: Fraction,
#                  perimeter: float,
#                  # border_triangles: Set[CachedPolygon],
#                  border_edges: Set[Segment[Fraction]]) -> None:
#         self.triangles = triangles
#         self.area = area
#         self.perimeter = perimeter
#         # self.border_triangles = border_triangles
#         self.border_edges = border_edges
#         self._hash = hash(tuple(triangles))
#
#         def __hash__(self) -> int:
#             return self._hash
#
#
# class Partition:
#     def __init__(self,
#                  triangular_map: Dict[CachedPolygon,
#                                       Dict[Segment[Fraction],
#                                            CachedPolygon]],
#                  chunk_map: Dict[Chunk, Dict[Set[Segment[Fraction]], Chunk]]
#                  ) -> None:
#         self.triangular_map = triangular_map
#         self.chunk_map = chunk_map
#         self._chunks = set(chunk_map)
#
#     @property
#     def chunks(self) -> Set[Chunk]:
#         return self._chunks
#
#     def is_articulation_chunk(self, chunk_to_remove: Chunk) -> bool:
#         """
#         Checks if removal of the given chunk
#         would disconnect the graph.
#         """
#         seen = set()
#         nextlevel = {next(chunk for chunk in self.chunk_map
#                           if chunk != chunk_to_remove)}
#         while nextlevel:
#             thislevel = nextlevel
#             nextlevel = set()
#             for chunk in thislevel:
#                 if chunk not in seen and chunk != chunk_to_remove:
#                     seen.add(chunk)
#                     nextlevel.update(self.chunk_map[chunk].values())
#         return len(self.chunk_map) - 1 != len(seen)
#
#     def is_articulation_triangle(self,
#                                  triangle_to_remove: CachedPolygon,
#                                  chunk_triangles: Set[CachedPolygon]
#                                  ) -> bool:
#         """
#         Checks if removal of the given triangle from the given set
#         would disconnect the corresponing graph.
#         """
#         seen = set()
#         nextlevel = {next(triangle for triangle in self.triangular_map
#                           if triangle != triangle_to_remove
#                           and triangle in chunk_triangles)}
#         while nextlevel:
#             thislevel = nextlevel
#             nextlevel = set()
#             for triangle in thislevel:
#                 if (triangle not in seen and triangle != triangle_to_remove
#                         and triangle in chunk_triangles):
#                     seen.add(triangle)
#                     nextlevel.update(self.triangular_map[triangle].values())
#         return len(chunk_triangles) - 1 != len(seen)
#
#     def move(self,
#              triangle: CachedPolygon,
#              *,
#              source: Chunk,
#              target: Chunk) -> None:
#         # assert triangle in source.border_triangles
#         source.triangles.remove(triangle)
#         target.triangles.add(triangle)
#         # source.border_triangles.remove(triangle)
#         # source.border_triangles.update(
#         #     neighbor for neighbor in self.triangular_map[triangle].values()
#         #     if neighbor in source)
#         # target.border_triangles.add(triangle)
#         # target.border_triangles -= set(
#         #     neighbor for neighbor in self.triangular_map[triangle].values()
#         #     if all(neighbors_neighbor in target.triangles
#         #            for neighbors_neighbor in self.triangular_map[neighbor].values()))
#         source.area -= triangle.area
#         target.area += triangle.area
#         source.perimeter += sum(
#             fast_length(edge) * (-1 if neighbor in source.triangles else 1)
#             for edge, neighbor in self.triangular_map[triangle].items())
#         target.perimeter += sum(
#             fast_length(edge) * (-1 if neighbor in target.triangles else 1)
#             for edge, neighbor in self.triangular_map[triangle].items())
#         source.border_edges -= {
#             edge for edge, neighbor in self.triangular_map[triangle].items()
#             if neighbor not in source.triangles}
#         source.border_edges |= {
#             edge for edge, neighbor in self.triangular_map[triangle].items()
#             if neighbor in source.triangles}
#         target.border_edges -= {
#             edge for edge, neighbor in self.triangular_map[triangle].items()
#             if neighbor in target.triangles}
#         target.border_edges |= {
#             edge for edge, neighbor in self.triangular_map[triangle].items()
#             if neighbor not in target.triangles}
#
#     def neighbors_of_chunk(self, chunk: Chunk) -> Iterator[Chunk]:
#         yield from self.chunk_map[chunk].values()
#
#     def neighbors_of_triangle(self, triangle: CachedPolygon
#                               ) -> Iterator[CachedPolygon]:
#         yield from self.triangular_map[triangle].values()
#
#     def edge_between_chunks(self, chunk: Chunk, other: Chunk
#                             ) -> FrozenSet[Segment[Fraction]]:
#         return next(edge for edge, neighbor in self.chunk_map[chunk].items()
#                     if neighbor == other)
#
#     def edge_between_triangles(self, triangle: CachedPolygon,
#                                other: CachedPolygon
#                                ) -> FrozenSet[Segment[Fraction]]:
#         return next(edge
#                     for edge, neighbor in self.triangular_map[triangle].items()
#                     if neighbor == other)
#
#     def unite(self, chunk: Chunk, other: Chunk) -> None:
#         assert chunk in self.chunk_map[other].values()
#         triangles = chunk.triangles | other.triangles
#         perimeter = (chunk.perimeter + other.perimeter
#                      - 2 * sum(map(fast_length, self.edge_between_chunks(
#                     chunk, other))))
#         # border_triangles_union = chunk.border_triangles | other.border_triangles
#         # border_triangles = border_triangles_union - {
#         #     triangle for triangle in border_triangles_union
#         #     if any(neighbor not in triangles
#         #            for neighbor in self.triangular_map[triangle].values())}
#         border_edges = chunk.border_edges ^ other.border_edges
#         union = Chunk(triangles=triangles,
#                       area=chunk.area + other.area,
#                       perimeter=perimeter,
#                       # border_triangles=border_triangles,
#                       border_edges=border_edges)
#         chunk_neighbors = {edge: neighbor
#                            for edge, neighbor in self.chunk_map[chunk].items()
#                            if neighbor != other}
#         other_neighbors = {edge: neighbor
#                            for edge, neighbor in self.chunk_map[other].items()
#                            if neighbor != chunk}
#         union_neighbors = {**chunk_neighbors, **other_neighbors}
#         self.chunk_map[union] = union_neighbors
#         for edge, neighbor in union_neighbors.items():
#             self.chunk_map[neighbor][edge] = union
#         del self.chunk_map[chunk]
#         del self.chunk_map[other]
#         self._chunks.remove(chunk)
#         self._chunks.remove(other)
#         self._chunks.add(union)
#
#
# # class Grid(Indexable[Coordinate], Shaped[Coordinate]):
# #     ___slots__ = '_elements',
# #
# #     def __init__(self, mapping: Dict[Shaped[Coordinate],
# #                                      Dict[Union[Segment[Coordinate],
# #                                                 List[Segment[Coordinate]]],
# #                                           Optional[Shaped[Coordinate]]]]
# #                  ) -> None:
# #         self._mapping = mapping
# #         self._elements = set(mapping)
# #         self._area = sum(element.area for element in mapping)
# #         self._border_elements = set(
# #             element for element, neighbor_per_edge in mapping.items()
# #             if None in neighbor_per_edge.values())
# #         self._border_polygons = set(
# #             polygon for element in self.border_elements
# #             for polygon in (element.border_polygons
# #                             if isinstance(element, Grid) else [element]))
# #         self._border_edges = set(
# #             edge for neighbor_per_edge in mapping.values()
# #             for edge, neighbor in neighbor_per_edge.items()
# #             if neighbor is None)
# #         self._perimeter = sum(map(fast_length, self._border_edges))
# #
# #     def __len__(self) -> int:
# #         return len(self._mapping)
# #
# #     @property
# #     def mapping(self):
# #         return self._mapping
# #
# #     @property
# #     def elements(self) -> Set[Shaped[Coordinate]]:
# #         return self._elements
# #
# #     @property
# #     def area(self) -> Coordinate:
# #         return self._area
# #
# #     @property
# #     def perimeter(self):
# #         return self._perimeter
# #
# #     @property
# #     def border_elements(self):
# #         return self._border_elements
# #
# #     @property
# #     def border_polygons(self):
# #         return self._border_polygons
# #
# #     def is_articulation_point(self,
# #                               element_to_remove: Shaped[Coordinate]) -> bool:
# #         """
# #         Checks if removal of the given element
# #         would disconnect the graph.
# #         """
# #         seen = set()
# #         nextlevel = {next(element for element in self._mapping
# #                           if element != element_to_remove)}
# #         while nextlevel:
# #             thislevel = nextlevel
# #             nextlevel = set()
# #             for element in thislevel:
# #                 if element not in seen and element != element_to_remove:
# #                     seen.add(element)
# #                     nextlevel.update(self._mapping[element])
# #         return len(self) - 1 == len(seen)
# #
# #     def remove(self, element: Shaped[Coordinate]) -> None:
# #         self._last_removed_element = element
# #         self._perimeter_before_removal = self._perimeter
# #         self._perimeter += sum(
# #             fast_length(edge) * (-1 if neighbor is None else 1)
# #             for edge, neighbor in self._mapping[element])
# #         self._border_polygons.remove(element)
# #         for neighbor in self._mapping[element].values():
# #             if neighbor not in self._border_polygons:
# #                 self._new_border_polygons.add(neighbor)
# #                 self._border_polygons.add(neighbor)
# #         self._last_removed_mapping = self._mapping[element]
# #         for edge, neighbor in self._mapping[element].items():
# #             self._mapping[neighbor][edge] = None
# #         del self._mapping[element]
# #         self._area -= element.area
# #
# #     def return_last_removed(self) -> None:
# #         self._perimeter = self._perimeter_before_removal
# #         for edge, neighbor in self._last_removed_mapping.items():
# #             self._mapping[neighbor][edge] = self._last_removed_element
# #         self._mapping[self._last_removed_element] = self._last_removed_mapping.copy()
# #         self.area += self._last_removed_element.area
# #         self._border_polygons.add(self._last_removed_element)
# #         self._border_polygons -= self._new_border_polygons
# #
# #     def neighbors_of(self, element: Shaped[Fraction]) -> Set[Shaped[Fraction]]:
# #         return set(self._mapping[element].values())
# #
# #     def edge_between(self,
# #                      element: Shaped[Fraction],
# #                      other: Shaped[Fraction]
# #                      ) -> Union[Segment[Fraction], List[Segment[Fraction]]]:
# #         return next(edge for edge, neighbor in self._mapping[element].items()
# #                     if neighbor == other)
#
#
# def to_partition(polygon: Polygon[Fraction],
#                  *,
#                  extra_points: Sequence[Point[Fraction]] = ()
#                  ) -> Partition:
#     """O(n²)"""
#     triangulation = Triangulation.constrained_delaunay(
#         polygon,
#         extra_points=extra_points,
#         context=context)
#     triangles = [CachedPolygon(contour, [])
#                  for contour in triangulation.triangles()]
#     triangles_per_edge = defaultdict(set)
#     chunks_per_edge = defaultdict(set)
#     for triangle in triangles:
#         chunk = Chunk(triangles={triangle},
#                       area=triangle.area,
#                       perimeter=fast_length(triangle.border),
#                       # border_triangles={triangle},
#                       border_edges=set(triangle.edges))
#         for edge in triangle.edges:
#             triangles_per_edge[edge].add(triangle)
#             chunks_per_edge[edge].add(chunk)
#     triangular_map = {}
#     chunk_map= {}
#     for edge, triangles in triangles_per_edge.items():
#         chunks = chunks_per_edge[edge]
#         for triangle, chunk in zip(triangles, chunks):
#             if triangle not in triangular_map:
#                 triangular_map[triangle] = {
#                     edge: next(iter(triangles - {triangle}), None)}
#             else:
#                 triangular_map[triangle][edge] = next(iter(triangles - {triangle}), None)
#             neighbor_chunk = next(iter(chunks - {chunk}), None)
#             if neighbor_chunk is not None:
#                 if chunk not in chunk_map:
#                     chunk_map[chunk] = {
#                         frozenset([edge]): neighbor_chunk}
#                 else:
#                     chunk_map[chunk][frozenset([edge])] = neighbor_chunk
#     return Partition(triangular_map=triangular_map,
#                      chunk_map=chunk_map)
#
#
# def chunk_split(polygon: Polygon[Fraction],
#                 *,
#                 requirement: Fraction,
#                 steiner_distance: float
#                 ) -> Tuple[Polygon[Fraction], Polygon[Fraction]]:
#     extra_points = steiner_points(polygon, delta=steiner_distance)
#     draw_polygon_with_extra_points(polygon, extra_points)
#     partition = to_partition(polygon, extra_points=extra_points)
#     draw_partition(partition); show()
#     draw_colored_partition(partition); show()
#     while True:
#         chunk = find_satisfactory_chunk(partition, requirement=requirement)
#         print(f"Chunk: {chunk.triangles if chunk is not None else None}")
#         if chunk is not None:
#             remainder = collect_remainder(chunk, partition)
#             if chunk.area > requirement:
#                 chunk, remainder = readjust(chunk,
#                                             remainder,
#                                             requirement=requirement,
#                                             partition=partition)
#             return unite(chunk.triangles), unite(remainder.triangles)
#         join_pairs(partition)
#         min_chunk_size = min(len(chunk.triangles)
#                              for chunk in partition.chunks)
#         if min_chunk_size >= 5:
#             swap_leaves(partition)
#
#
# def find_satisfactory_chunk(partition: Partition,
#                             *,
#                             requirement: Fraction) -> Optional[Chunk]:
#     large_chunks = (chunk for chunk in partition.chunks
#                     if chunk.area >= requirement
#                     and not partition.is_articulation_chunk(chunk))
#     return min(large_chunks, key=attrgetter('perimeter'), default=None)
#
#
# def readjust(chunk: Chunk,
#              neighbor: Chunk,
#              *,
#              requirement: Fraction,
#              partition: Partition) -> Tuple[Chunk, Chunk]:
#     draw_partition(partition); draw_chunk(chunk, color='blue', fill=True); show()
#     while True:
#         border_triangles = {
#             triangle for triangle in chunk.triangles
#             if any(neighbor_triangle in neighbor.triangles
#                    for neighbor_triangle in partition.triangular_map[triangle].values())
#             and not partition.is_articulation_triangle(triangle, chunk.triangles)}
#         if not border_triangles:
#             raise RuntimeError
#         while border_triangles:
#             triangle = border_triangles.pop()
#             partition.move(triangle, source=chunk, target=neighbor)
#             draw_partition(partition); draw_chunk(chunk, color='blue', fill=True); show()
#             if chunk.area == requirement:
#                 return chunk, neighbor
#             elif chunk.area < requirement:
#                 partition.move(triangle, source=neighbor, target=chunk)
#                 node_to_split = triangle
#                 pivot_index = find_pivot(node_to_split, partition=partition, chunk=chunk)
#                 vertices = node_to_split.border.vertices
#                 if DRAW:
#                     draw_partition(partition); draw_chunk(chunk, color='blue', fill=True); draw_triangle(node_to_split, color='red', fill=True);
#                     draw(vertices[pivot_index]); show()
#                 new_vertex = inverse_shoelace(requirement + triangle.area - chunk.area,
#                                               endpoint=vertices[pivot_index],
#                                               base_start=vertices[pivot_index - 1],
#                                               base_end=vertices[pivot_index - 2])
#                 if DRAW:
#                     draw_partition(partition); draw_chunk(chunk, color='blue', fill=True); draw_triangle(node_to_split, color='red', fill=True);
#                     draw(vertices[pivot_index]); draw(new_vertex); show()
#                 divide(partition, node=node_to_split,
#                        pivot_index=pivot_index,
#                        new_vertex=new_vertex,
#                        chunk=chunk,
#                        neighbor=neighbor)
#                 return chunk, neighbor
#
#
# def collect_remainder(chunk: Chunk, partition: Partition) -> Chunk:
#     while len(partition.chunks) > 2:
#         neighbor = next(partition.neighbors_of_chunk(chunk))
#         neighbors_neighbor = next(
#             node for node in partition.neighbors_of_chunk(neighbor)
#             if node != chunk)
#         partition.unite(neighbor, neighbors_neighbor)
#     return next(partition.neighbors_of_chunk(chunk))
#
#
# def join_pairs(partition: Partition) -> None:
#     smallest_chunk_size = min(len(chunk.triangles)
#                               for chunk in partition.chunks)
#     print(f"Smallest chunk size: {smallest_chunk_size}")
#     potential_sources = {chunk for chunk in partition.chunks
#                          if len(chunk.triangles) == smallest_chunk_size}
#     potential_targets = set(partition.chunks)
#     while potential_sources:
#         merge_least_compact(partition,
#                             potential_sources,
#                             potential_targets)
#         print(len(potential_sources), len(potential_targets))
#
#
# def find_leaf_pair(graph: nx.Graph
#                    ) -> Optional[Tuple[FrozenSet[Polygon[Fraction]],
#                                        FrozenSet[Polygon[Fraction]]]]:
#     for node in graph:
#         neighbors = list(nx.neighbors(graph, node))
#         if len(neighbors) == 1:
#             neighbor = neighbors[0]
#             if len(list(nx.neighbors(graph, neighbor))) == 2:
#                 return node, neighbor
#     return None
#
#
# def merge_least_compact(
#         partition: Partition,
#         potential_sources: Set[Chunk],
#         potential_targets: Set[Chunk]) -> None:
#     source = min(potential_sources,
#                  key=lambda chunk: chunk.area ** 0.5 / chunk.perimeter)
#     # draw_chunks_graph(graph); draw_source_node(source); show()
#     best_target = None
#     max_compactness = 0
#     neighbors = [neighbor for neighbor in partition.neighbors_of_chunk(source)
#                  if neighbor in potential_targets]
#     if not neighbors:
#         potential_sources.remove(source)
#         return
#     for neighbor in neighbors:
#         area = source.area + neighbor.area
#         perimeter = (source.perimeter + neighbor.perimeter
#                      - 2 * sum(map(fast_length, partition.edge_between_chunks(
#                     source, neighbor))))
#         compactness = area ** 0.5 / perimeter
#         if compactness > max_compactness:
#             best_target = neighbor
#             max_compactness = compactness
#     potential_sources -= {source, best_target}
#     potential_targets -= {source, best_target}
#     if len(source.triangles) > PLOT_STARTING_NUMBER:
#         draw_colored_partition(partition); draw_source_chunk(source); draw_target_chunk(best_target); show()
#     partition.unite(source, best_target)
#
#
# def union_attributes(chunk: FrozenSet[Polygon[Fraction]],
#                      other: FrozenSet[Polygon[Fraction]],
#                      *,
#                      graph: nx.Graph) -> float:
#     area = graph.nodes[chunk]['area'] + graph.nodes[other]['area']
#     common_multisegment = graph[chunk][other]['edge']
#     perimeter = (graph.nodes[chunk]['perimeter']
#                  + graph.nodes[other]['perimeter']
#                  - 2 * fast_length(common_multisegment))
#     compactness = float(area) ** 0.5 / perimeter
#     return area, perimeter, compactness
#
#
# def swap_leaves(partition: Partition) -> None:
#     for chunk in partition.chunks:
#         # for triangle in chunk.border_triangles:
#         for triangle in set(chunk.triangles):
#             triangle_neighbors = set(
#                 neighbor
#                 for neighbor in partition.neighbors_of_triangle(triangle)
#                 if neighbor is not None)
#             if len(triangle_neighbors) != 3:
#                 continue
#             if sum(neighbor in chunk.triangles
#                    for neighbor in triangle_neighbors) != 1:
#                 continue
#             chunk_triangle = next(neighbor for neighbor in triangle_neighbors
#                                   if neighbor in chunk.triangles)
#             first_neighbor, second_neighbor = triangle_neighbors - {chunk_triangle}
#             chunk_neighbors = set(partition.neighbors_of_chunk(chunk))
#             first_neighbor_chunk = next(
#                 chunk_neighbor for chunk_neighbor in chunk_neighbors
#                 if first_neighbor in chunk_neighbor.triangles)
#             second_neighbor_chunk = next(
#                 chunk_neighbor for chunk_neighbor in chunk_neighbors
#                 if second_neighbor in chunk_neighbor.triangles)
#             if first_neighbor_chunk != second_neighbor_chunk:
#                 continue
#             partition.move(triangle,
#                            source=chunk,
#                            target=first_neighbor_chunk)
#             break
#
#
# def inverse_shoelace(area: Fraction,
#                      endpoint: Point,
#                      base_start: Point,
#                      base_end: Point) -> Point:
#     orientation = Contour([endpoint, base_start, base_end]).orientation
#     sign = 1 if orientation is Orientation.COUNTERCLOCKWISE else -1
#     if base_end.x - base_start.x == 0:
#         numerator = 2 * area + sign * (base_start.x * base_start.y
#                                        - endpoint.x * base_start.y)
#         denominator = sign * (base_start.x - endpoint.x)
#         y = numerator / denominator
#         return Point(base_start.x, y)
#     slope, intercept = slope_intercept(base_start, base_end)
#     x_diff = base_start.x - endpoint.x
#     numerator = 2 * area + sign * (base_start.x * endpoint.y
#                                    - endpoint.x * base_start.y
#                                    - intercept * x_diff)
#     denominator = sign * (endpoint.y - base_start.y + slope * x_diff)
#     x = numerator / denominator
#     y = slope * x + intercept
#     return Point(x, y)
#
#
# def slope_intercept(start: Point, end: Point) -> Tuple[Fraction, Fraction]:
#     slope = (end.y - start.y) / (end.x - start.x)
#     intercept = end.y - slope * end.x
#     return slope, intercept
#
#
# def compactness(*,
#                 area_: float,
#                 perimeter: float) -> float:
#     return area_ ** 0.5 / perimeter
#
#
# def unite(geometries: Iterable[Compound]) -> Compound:
#     return reduce(or_, geometries, EMPTY)
#
#
# def find_pivot(triangle: CachedPolygon,
#                *,
#                partition: Partition,
#                chunk: Chunk) -> int:
#     predecessor = next(
#         (neighbor for neighbor in partition.neighbors_of_triangle(triangle)
#          if neighbor in chunk.triangles), None)
#     if predecessor is not None:
#         unoriented_edge = partition.edge_between_triangles(triangle,
#                                                            predecessor)
#         oriented_edge = next(edge for edge in triangle.edges
#                              if edge == unoriented_edge)
#         return triangle.border.vertices.index(oriented_edge.end)
#     neighbor = next(
#         (neighbor for neighbor in partition.neighbors_of_triangle(triangle)
#          if neighbor not in chunk.triangles), None)
#     if neighbor is None:
#         return 0  # any index will do
#     unoriented_edge = partition.edge_between_triangles(triangle, neighbor)
#     oriented_edge = next(edge for edge in triangle.edges
#                          if edge == unoriented_edge)
#     return triangle.border.vertices.index(oriented_edge.start)
#
#
# def divide(partition: Partition,
#            node: CachedPolygon,
#            pivot_index: int,
#            new_vertex: Point,
#            chunk: Chunk,
#            neighbor: Chunk) -> None:
#     vertices = node.border.vertices
#     pivot_vertex = vertices[pivot_index]
#     left_vertex = vertices[pivot_index - 1]
#     right_vertex = vertices[pivot_index - 2]
#     left_segment = Segment(pivot_vertex, left_vertex)
#     right_segment = Segment(right_vertex, pivot_vertex)
#     third_segment = Segment(left_vertex, right_vertex)
#     left_part = CachedPolygon(Contour([pivot_vertex,
#                                        left_vertex,
#                                        new_vertex]), [])
#     right_part = CachedPolygon(Contour([right_vertex,
#                                         pivot_vertex,
#                                         new_vertex]), [])
#     left_neighbor = partition.triangular_map[node][left_segment]
#     right_neighbor = partition.triangular_map[node][right_segment]
#     third_neighbor = partition.triangular_map[node][third_segment]
#     partition.triangular_map[left_part] = {Segment(pivot_vertex, new_vertex): right_part}
#     partition.triangular_map[right_part] = {Segment(pivot_vertex, new_vertex): left_part}
#     if left_neighbor is not None:
#         partition.triangular_map[left_part][left_segment] = left_neighbor
#         partition.triangular_map[left_neighbor][left_segment] = left_part
#     if right_neighbor is not None:
#         partition.triangular_map[right_part][right_segment] = right_neighbor
#         partition.triangular_map[right_neighbor][right_segment] = right_part
#     if third_neighbor is not None:
#         third_neighbor_only_vertex = next(
#             vertex for vertex in third_neighbor.border.vertices
#             if vertex not in {left_vertex, right_vertex})
#         third_neighbor_left_part = CachedPolygon(
#             Contour([new_vertex, left_vertex, third_neighbor_only_vertex]), [])
#         third_neighbor_right_part = CachedPolygon(
#             Contour([right_vertex, new_vertex, third_neighbor_only_vertex]),
#             [])
#         partition.triangular_map[third_neighbor_left_part] = {Segment(new_vertex, third_neighbor_only_vertex): third_neighbor_right_part}
#         partition.triangular_map[third_neighbor_right_part] = {Segment(new_vertex, third_neighbor_only_vertex): third_neighbor_left_part}
#         partition.triangular_map[left_part][Segment(left_vertex, new_vertex)] = third_neighbor_left_part
#         partition.triangular_map[third_neighbor_left_part][Segment(left_vertex, new_vertex)] = left_part
#         partition.triangular_map[right_part][Segment(right_vertex, new_vertex)] = third_neighbor_right_part
#         partition.triangular_map[third_neighbor_right_part][Segment(right_vertex, new_vertex)] = right_part
#         third_neighbor_left_segment = Segment(left_vertex,
#                                               third_neighbor_only_vertex)
#         third_neighbor_right_segment = Segment(right_vertex,
#                                                third_neighbor_only_vertex)
#         third_neighbor_left_neighbor = next(
#             (neighbor for neighbor in partition.neighbors_of_triangle(third_neighbor)
#              if partition.edge_between_triangles(third_neighbor, neighbor) == third_neighbor_left_segment), None)
#         third_neighbor_right_neighbor = next(
#             (neighbor for neighbor in partition.neighbors_of_triangle(third_neighbor)
#              if partition.edge_between_triangles(third_neighbor, neighbor) == third_neighbor_right_segment), None)
#         if third_neighbor_left_neighbor is not None:
#             partition.triangular_map[third_neighbor_left_part][third_neighbor_left_segment] = third_neighbor_left_neighbor
#             partition.triangular_map[third_neighbor_left_neighbor][third_neighbor_left_segment] = third_neighbor_left_part
#         if third_neighbor_right_neighbor is not None:
#             partition.triangular_map[third_neighbor_right_part][third_neighbor_right_segment] = third_neighbor_right_neighbor
#             partition.triangular_map[third_neighbor_right_neighbor][third_neighbor_right_segment] = third_neighbor_right_part
#         del partition.triangular_map[third_neighbor]
#         chunk_with_third_neighbor = next(chunk_ for chunk_ in (chunk, neighbor)
#                                          if third_neighbor in chunk_.triangles)
#         chunk_with_third_neighbor.triangles.remove(third_neighbor)
#         chunk_with_third_neighbor.triangles.update((third_neighbor_left_part,
#                                                     third_neighbor_right_part))
#
#     del partition.triangular_map[node]
#     chunk.triangles.remove(node)
#     chunk.triangles.add(left_part)
#     neighbor.triangles.add(right_part)
#
#
# def draw_progress(accumulated_nodes: Set[CachedPolygon],
#                   current_polygon: CachedPolygon) -> None:
#     if DRAW:
#         for node in accumulated_nodes:
#             draw(node, fill=True, color='green', alpha=0.3)
#         draw(current_polygon, fill=True, color='green')
#
#
# def draw_border(border_nodes: Set[CachedPolygon]) -> None:
#     if DRAW:
#         for node in border_nodes:
#             draw(node, fill=True, color='yellow', alpha=0.7)
#
#
# def add_title_and_show(title: str) -> None:
#     if DRAW:
#         plt.title(title)
#         show()
#
#
# def draw_split(graph: nx.Graph,
#                accumulated_nodes: Set[CachedPolygon],
#                parts: Tuple[CachedPolygon, CachedPolygon]) -> None:
#     if DRAW:
#         for node in accumulated_nodes:
#             draw(node, fill=True, color='green', alpha=0.3)
#         draw(parts[0], fill=True, color='green')
#         draw(parts[1], fill=True, color='orange')
#         for node in graph:
#             draw(node, color='black')
#         show()
#
#
# def fast_length(geometry: Linear[Fraction]) -> float:
#     original_sqrt = context.sqrt
#     context._sqrt = math.sqrt
#     result = geometry.length
#     context._sqrt = original_sqrt
#     return result
#
#
# def fast_area(nodes: Iterable[Polygon]) -> float:
#     original_sqrt = context.sqrt
#     context._sqrt = math.sqrt
#     result = sum(node.area for node in nodes)
#     context._sqrt = original_sqrt
#     return result
#
#
# def steiner_points(polygon: Polygon[Fraction],
#                    delta: Fraction) -> List[Point[Fraction]]:
#     xs = [vertex.x for vertex in polygon.border.vertices]
#     ys = [vertex.y for vertex in polygon.border.vertices]
#     x_start = min(xs)
#     y_start = min(ys)
#     width = max(xs) - x_start
#     height = max(ys) - y_start
#     xs_grid_count = width // delta
#     ys_grid_count = height // delta
#     grid = [Point(x_start + x * delta, y_start + y * delta)
#             for x, y in product(range(1, xs_grid_count),
#                                 range(1, ys_grid_count))]
#     constrained_grid = [point for point in grid if point in polygon]
#     points_on_border = []
#     for segment in polygon.border.segments:
#         extra_points_count = Fraction(fast_length(segment)) // delta
#         if not extra_points_count:
#             continue
#         dx = (segment.end.x - segment.start.x) / extra_points_count
#         dy = (segment.end.y - segment.start.y) / extra_points_count
#         for i in range(1, extra_points_count + 1):
#             points_on_border.append(Point(segment.start.x + i * dx,
#                                           segment.start.y + i * dy))
#     return constrained_grid + points_on_border
#
#
# def centroid(node: Polygon) -> Point:
#     original_sqrt = context.sqrt
#     context._sqrt = math.sqrt
#     result = node.centroid
#     context._sqrt = original_sqrt
#     return result
#
#
# def draw_polygon_with_extra_points(polygon: Polygon[Fraction],
#                                    extra_points: Iterable[Point[Fraction]]
#                                    ) -> None:
#     if DRAW:
#         draw(polygon)
#         for point in extra_points:
#             draw(point, color='b')
#         show()
#
#
# def draw_partition(partition: Partition,
#                    *,
#                    color: str = 'black',
#                    fill: bool = False) -> None:
#     if DRAW:
#         for chunk in partition.chunks:
#             draw_chunk(chunk, color=color, fill=fill)
#
#
# def draw_chunk(chunk: Chunk,
#                *,
#                color: str = 'black',
#                fill: bool = False) -> None:
#     if DRAW:
#         for triangle in chunk.triangles:
#             draw(triangle, color=color, fill=fill)
#
#
# def draw_triangle(triangle: CachedPolygon,
#                   *,
#                   color: str = 'black',
#                   fill: bool = False) -> None:
#     if DRAW:
#         draw(triangle, color=color, fill=fill)
#
#
# def draw_colored_partition(partition: Partition) -> None:
#     if DRAW:
#         for chunk in partition.chunks:
#             color = (random.uniform(0, 1),
#                      random.uniform(0, 1),
#                      random.uniform(0, 1))
#             for triangle in chunk.triangles:
#                 draw(triangle, fill=True, color=color)
#                 draw(triangle, color='k')
#
#
# def draw_source_chunk(chunk: Chunk) -> None:
#     if DRAW:
#         for triangle in chunk.triangles:
#             draw(triangle, color='red', linewidth=2)
#
#
# def draw_target_chunk(chunk: Chunk) -> None:
#     if DRAW:
#         for triangle in chunk.triangles:
#             draw(triangle, color='blue', linewidth=2)
