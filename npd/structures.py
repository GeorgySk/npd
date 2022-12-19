"""Various classes designed to simplify the algorithm"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import chain
from typing import (Dict,
                    FrozenSet,
                    Iterator,
                    List,
                    Optional,
                    Sequence)

from gon.base import (Contour,
                      Point,
                      Polygon,
                      Segment,
                      Triangulation)
from gon.hints import Scalar
from ground.base import get_context

from npd.utils import fast_length

NeighborChunkPerEdge = Dict[FrozenSet[Segment[Fraction]], 'Chunk']
TriangularMap = Dict['CachedPolygon', Dict[Segment[Fraction], 'CachedPolygon']]


@dataclass
class ConvexPartition:
    """
    Class to keep information on domains, their corresponding
    countersegments, vertices on each side, and the area difference
    """
    domain: Segment[Fraction]
    countersegment: Segment[Fraction]
    right_vertices: List[Point[Fraction]]
    left_vertices: [Point[Fraction]]
    area_difference: Fraction


class CachedPolygon(Polygon):
    """
    Polygon with hashed perimeter
    """
    slots = '_hash', '_perimeter'

    def __init__(self,
                 border: Contour[Scalar],
                 holes: Optional[Sequence[Contour[Scalar]]] = None) -> None:
        if holes is None:
            holes = []
        super().__init__(border, holes)
        self._hash = super().__hash__()
        self._perimeter = (fast_length(self.border)
                           + sum(map(fast_length, holes)))

    def __hash__(self) -> int:
        return self._hash

    @property
    def perimeter(self) -> float:
        """
        :return: hashed perimeter
        """
        return self._perimeter


class Chunk:
    """
    Stores collection of triangles with their total area and perimeter
    """
    def __init__(self, *,
                 triangles: Dict[CachedPolygon, None],
                 area: Fraction,
                 perimeter: float) -> None:
        self.triangles = triangles
        self.area = area
        self._perimeter = perimeter
        if not isinstance(perimeter, float):
            raise TypeError("Bad perimeter value")
        self._hash = hash(tuple(triangles))

    def __hash__(self) -> int:
        return self._hash

    @property
    def perimeter(self) -> float:
        return self._perimeter

    @perimeter.setter
    def perimeter(self, value) -> None:
        if not isinstance(value, float):
            raise TypeError("Got unexpected type")
        self._perimeter = value

    def compactness(self) -> float:
        """
        :return: compactness of the union of triangles
        """
        return self.area ** 0.5 / self.perimeter


class Partition:
    """
    Class to keep relation between different chunks and triangles
    the polygon is split into.
    """
    def __init__(self,
                 polygon: Polygon,
                 *,
                 extra_points: Sequence[Point[Fraction]] = ()) -> None:
        """O(n²)"""
        triangulation = Triangulation.constrained_delaunay(
            polygon,
            extra_points=extra_points,
            context=get_context())
        triangles = list(map(CachedPolygon, triangulation.triangles()))
        self.triangular_map = to_triangular_map(triangles)
        self.chunk_map = to_chunk_map(self.triangular_map)
        self.total_area = polygon.area

    def is_articulation_chunk(self, chunk_to_remove: Chunk) -> bool:
        """
        Checks if removal of the given chunk
        would disconnect the graph.
        """
        root = next(chunk for chunk in self.chunk_map
                    if chunk != chunk_to_remove)
        seen, stack = {chunk_to_remove}, [root]
        while stack:
            chunk = stack.pop()
            seen.add(chunk)
            stack.extend(neighbor
                         for neighbor in self.chunk_map[chunk].values()
                         if neighbor not in seen)
        return len(self.chunk_map) != len(seen)

    def is_articulation_triangle(self,
                                 triangle_to_remove: CachedPolygon,
                                 chunk_triangles: Dict[CachedPolygon, None]
                                 ) -> bool:
        """
        Checks if removal of the given triangle from the given set
        would disconnect the corresponding graph.
        """
        if len(chunk_triangles) == 1:
            return False
        root = next(triangle for triangle in self.triangular_map
                    if triangle != triangle_to_remove
                    and triangle in chunk_triangles)
        seen, stack = {triangle_to_remove}, [root]
        while stack:
            chunk = stack.pop()
            seen.add(chunk)
            stack.extend(neighbor
                         for neighbor in self.triangular_map[chunk].values()
                         if neighbor not in seen
                         and neighbor in chunk_triangles)
        return len(chunk_triangles) != len(seen)

    def move(self,
             triangle: CachedPolygon,
             *,
             source: Chunk,
             target: Chunk) -> None:
        """
        Moves the given triangle from the source chunk to the target
        chunk
        :param triangle: will be moved from one chunk to another
        :param source: chunk that will lose the triangle
        :param target: chunk that will receive the triangle
        """
        source.triangles.pop(triangle)
        target.triangles.setdefault(triangle)
        source.area -= triangle.area
        target.area += triangle.area
        source.perimeter += sum(
            fast_length(edge) * (1 if neighbor in source.triangles else -1)
            for edge, neighbor in self.triangular_map[triangle].items())
        if not isinstance(source.perimeter, float):
            raise TypeError("Something went wrong when changing perimeter")
        target.perimeter += sum(
            fast_length(edge) * (-1 if neighbor in target.triangles else 1)
            for edge, neighbor in self.triangular_map[triangle].items())
        if not isinstance(target.perimeter, float):
            raise TypeError("Something went wrong when changing perimeter")

    def chunk_neighbors(self, chunk: Chunk) -> List[Chunk]:
        """
        :param chunk: any chunk in the partition
        :return: a list of neighbor chunks
        """
        return list(self.chunk_map[chunk].values())

    def triangle_neighbors(self, triangle: CachedPolygon
                           ) -> List[CachedPolygon]:
        """
        :param triangle: any triangle in the triangular map
        :return: a list of triangle's neighbors
        """
        return list(neighbor
                    for neighbor in self.triangular_map[triangle].values()
                    if neighbor is not None)

    def edge_between_chunks(self,
                            chunk: Chunk,
                            other: Chunk
                            ) -> FrozenSet[Segment[Fraction]]:
        """
        :param chunk: chunk in the partition
        :param other: neighbor chunk
        :return: edge between the given chunks as a set of segments
        """
        return next(edge for edge, neighbor in self.chunk_map[chunk].items()
                    if neighbor == other)

    def edge_between_triangles(self,
                               triangle: CachedPolygon,
                               other: CachedPolygon
                               ) -> Segment[Fraction]:
        """
        :param triangle: triangle in the triangular map
        :param other: neighbor of the triangle
        :return: segment that the given triangles share
        """
        return next(edge
                    for edge, neighbor in self.triangular_map[triangle].items()
                    if neighbor == other)

    def unite(self, chunk: Chunk, other: Chunk) -> Chunk:
        """
        Creates a new chunk out of two given ones, replaces them with
        the new one, and updates all the relations.
        :param chunk: any chunk from the partition
        :param other: neighbor chunk
        :return: newly created chunk as their union
        """
        union = Chunk(triangles={**chunk.triangles, **other.triangles},
                      area=chunk.area + other.area,
                      perimeter=self.union_perimeter(chunk, other))
        self.chunk_map[union] = self._union_neighbors_per_edge(chunk, other)
        self._add_chunk_to_neighbors(union)
        self._remove_chunk_from_neighbors(chunk)
        self._remove_chunk_from_neighbors(other)
        self.chunk_map.pop(chunk)
        self.chunk_map.pop(other)
        return union

    def union_perimeter(self, chunk: Chunk, other: Chunk) -> float:
        """
        :param chunk: any chunk from the partition
        :param other: neighbor chunk
        :return: perimeter of a chunk that could be obtained if the
        chunks were to be united
        """
        edge = self.edge_between_chunks(chunk, other)
        edge_length = sum(map(fast_length, edge))
        return chunk.perimeter + other.perimeter - 2 * edge_length

    def _union_neighbors_per_edge(self,
                                  chunk: Chunk,
                                  other: Chunk) -> NeighborChunkPerEdge:
        chunk_neighbors = {neighbor: edge
                           for edge, neighbor in self.chunk_map[chunk].items()
                           if neighbor != other}
        other_neighbors = {neighbor: edge
                           for edge, neighbor in self.chunk_map[other].items()
                           if neighbor != chunk}
        neighbors = dict.fromkeys(chain(chunk_neighbors, other_neighbors))
        return {chunk_neighbors.get(neighbor, frozenset())
                | other_neighbors.get(neighbor, frozenset()): neighbor
                for neighbor in neighbors}

    def _add_chunk_to_neighbors(self, chunk: Chunk) -> None:
        if chunk not in self.chunk_map:
            raise KeyError("Neighbors-per-edge mapping of the given chunk "
                           "must be added to chunk map before trying to add "
                           "the chunk in its neighbors mappings")
        for edge, neighbor in self.chunk_map[chunk].items():
            self.chunk_map[neighbor][edge] = chunk

    def _remove_chunk_from_neighbors(self, chunk: Chunk) -> None:
        for neighbor in self.chunk_neighbors(chunk):
            self.chunk_map[neighbor] = {
                edge: node for edge, node in self.chunk_map[neighbor].items()
                if node != chunk}

    def is_touching_other_chunk(self,
                                triangle: CachedPolygon,
                                *,
                                this_chunk: Chunk) -> bool:
        """
        Checks if the given triangle shares a border with any chunk
        other than the one this triangle belongs to
        :param triangle: any triangle from the triangular map
        :param this_chunk: chunk this triangle belongs to
        :return: True if the condition satisfies, False otherwise
        """
        return any(neighbor not in this_chunk.triangles
                   for neighbor in self.triangle_neighbors(triangle))

    def border_triangles(self, chunk: Chunk) -> Iterator[CachedPolygon]:
        """
        :param chunk: any chunk from the partition
        :return: iterator over triangles on the border of the chunk
        """
        for triangle in chunk.triangles:
            if self.is_touching_other_chunk(triangle, this_chunk=chunk):
                yield triangle

    def movable_triangles(self, chunk: Chunk) -> Iterator[CachedPolygon]:
        """
        :return: iterator over triangles of the given chunk that can be
        moved without causing discontinuity
        """
        yield from (triangle for triangle in self.border_triangles(chunk)
                    if not self.is_articulation_triangle(triangle,
                                                         chunk.triangles))


def to_triangular_map(triangles: List[CachedPolygon]) -> TriangularMap:
    """
    :param triangles: a list of triangles
    :return: a map that stores relation between neighboring triangles
    and shared segments
    """
    triangular_map = {}
    triangle_by_edge = {}
    for triangle in triangles:
        triangular_map[triangle] = {}
        for edge in triangle.edges:
            neighbor = triangle_by_edge.get(edge)
            triangular_map[triangle][edge] = neighbor
            if neighbor is None:
                triangle_by_edge[edge] = triangle
            else:
                triangular_map[neighbor][edge] = triangle
    return triangular_map


def to_chunk_map(triangular_map: TriangularMap) -> Dict[Chunk,
                                                        NeighborChunkPerEdge]:
    """
    :param triangular_map: input triangular map
    :return: same triangular map but where triangles are wrapped as
    chunks and segments wrapped as sets of single segment
    """
    chunk_per_triangle = {triangle: to_chunk(triangle)
                          for triangle in triangular_map}
    return {chunk_per_triangle[triangle]:
            {frozenset([edge]): chunk_per_triangle[neighbor]
             for edge, neighbor in neighbors_map.items()
             if neighbor is not None}
            for triangle, neighbors_map in triangular_map.items()}


def to_chunk(triangle: CachedPolygon) -> Chunk:
    """
    Wraps triangle as a chunk.
    :param triangle: input triangle
    :return: chunk wrapping the triangle with information about area
    and the perimeter
    """
    return Chunk(triangles={triangle: None},
                 area=triangle.area,
                 perimeter=fast_length(triangle.border))
