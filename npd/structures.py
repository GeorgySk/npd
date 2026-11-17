"""Various classes designed to simplify the algorithm."""
from collections.abc import (Collection,
                             Iterator,
                             Sequence)
from dataclasses import dataclass
from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon,
                      Segment,
                      Triangulation)
from gon.hints import Scalar
from ground.base import get_context

from npd.utils import fast_length

SharedEdgesPerNeighbor = dict["Chunk", frozenset[Segment[Fraction]]]
TriangularMap = dict["CachedPolygon",
                     dict[Segment[Fraction], "CachedPolygon | None"]]


@dataclass
class ConvexPartition:
    """Domain with its countersegment, side vertices and area difference."""  # noqa: D204
    domain: Segment[Fraction]
    countersegment: Segment[Fraction]
    right_vertices: list[Point[Fraction]]
    left_vertices: list[Point[Fraction]]
    area_difference: Fraction


class CachedPolygon(Polygon):  # type: ignore[misc]
    """Polygon with hashed perimeter."""

    __slots__ = "_hash", "_perimeter"

    def __init__(self,
                 border: Contour[Scalar],
                 holes: Sequence[Contour[Scalar]] | None = None) -> None:
        """Creates a polygon, computing its hash and perimeter only once."""
        if holes is None:
            holes = []
        super().__init__(border, holes)
        self._hash: int = super().__hash__()
        self._perimeter = (fast_length(self.border)
                           + sum(map(fast_length, holes)))

    def __hash__(self) -> int:
        """Hash computed at creation instead of on every call."""
        return self._hash

    @property
    def perimeter(self) -> float:
        """Hashed perimeter of the polygon."""
        return self._perimeter


@dataclass(eq=False, slots=True)
class Chunk:
    """
    Triangles with their total area and the perimeter of their union.

    Area and perimeter are cached and updated incrementally as triangles
    move between chunks. ``eq=False`` keeps identity equality and hashing,
    so a chunk stays a valid dict key while its triangles change.
    """

    triangles: dict[CachedPolygon, None]
    area: Fraction
    perimeter: float

    def compactness(self) -> float:
        """Compactness of the union of triangles."""
        return self.area ** 0.5 / self.perimeter


class Partition:
    """Relations between the chunks and triangles a polygon is split into."""

    def __init__(self,
                 polygon: Polygon,
                 *,
                 extra_points: Sequence[Point[Fraction]] = ()) -> None:
        """O(n²)."""
        triangulation = Triangulation.constrained_delaunay(
            polygon,
            extra_points=extra_points,
            context=get_context())
        triangles = list(map(CachedPolygon, triangulation.triangles()))
        self.triangular_map = to_triangular_map(triangles)
        self.chunk_map = to_chunk_map(self.triangular_map)
        self.total_area = polygon.area

    def is_articulation_chunk(self, chunk_to_remove: Chunk) -> bool:
        """Checks if removal of the given chunk would disconnect the graph."""
        root = next(chunk for chunk in self.chunk_map
                    if chunk != chunk_to_remove)
        seen, stack = {chunk_to_remove}, [root]
        while stack:
            chunk = stack.pop()
            seen.add(chunk)
            stack.extend(neighbor
                         for neighbor in self.chunk_map[chunk]
                         if neighbor not in seen)
        return len(self.chunk_map) != len(seen)

    def is_articulation_triangle(self,
                                 triangle_to_remove: CachedPolygon,
                                 chunk_triangles: Collection[CachedPolygon]
                                 ) -> bool:
        """Checks if removing the triangle disconnects the given triangles."""
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
                         if neighbor is not None
                         and neighbor not in seen
                         and neighbor in chunk_triangles)
        return len(chunk_triangles) != len(seen)

    def move(self,
             triangle: CachedPolygon,
             *,
             source: Chunk,
             target: Chunk) -> None:
        """
        Moves the given triangle from the source chunk to the target chunk.

        Also updates the edges the chunks share in the chunk map: each
        edge of the triangle stops being shared between the source and the
        chunk on its other side, and starts being shared between the
        target and that chunk.
        :param triangle: will be moved from one chunk to another
        :param source: chunk that will lose the triangle
        :param target: chunk that will receive the triangle
        """
        # Any chunk on the other side of the triangle's edges is the
        # source itself or one of its neighbors (the target among them).
        candidates = [source, *self.chunk_neighbors(source)]
        for edge, neighbor in self.triangular_map[triangle].items():
            if neighbor is None:
                continue
            owner = next(chunk for chunk in candidates
                         if neighbor in chunk.triangles)
            if owner is not source:
                self._remove_shared_edge(source, owner, edge)
            if owner is not target:
                self._add_shared_edge(target, owner, edge)
        source.triangles.pop(triangle)
        target.triangles.setdefault(triangle)
        source.area -= triangle.area
        target.area += triangle.area
        source.perimeter += sum(
            fast_length(edge) * (1 if neighbor in source.triangles else -1)
            for edge, neighbor in self.triangular_map[triangle].items())
        target.perimeter += sum(
            fast_length(edge) * (-1 if neighbor in target.triangles else 1)
            for edge, neighbor in self.triangular_map[triangle].items())

    def chunk_neighbors(self, chunk: Chunk) -> list[Chunk]:
        """
        Lists the neighbor chunks of the given chunk.

        :param chunk: any chunk in the partition
        :return: a list of neighbor chunks
        """
        return list(self.chunk_map[chunk])

    def triangle_neighbors(self, triangle: CachedPolygon
                           ) -> list[CachedPolygon]:
        """
        Lists the neighbor triangles of the given triangle.

        :param triangle: any triangle in the triangular map
        :return: a list of triangle's neighbors
        """
        return [neighbor
                for neighbor in self.triangular_map[triangle].values()
                if neighbor is not None]

    def shared_edges(self,
                     chunk: Chunk,
                     other: Chunk) -> frozenset[Segment[Fraction]]:
        """
        Finds the segments shared by two chunks.

        :param chunk: chunk in the partition
        :param other: chunk in the partition
        :return: segments on the border between the given chunks, empty
        if they aren't neighbors
        """
        return self.chunk_map[chunk].get(other, frozenset())

    def edge_between_triangles(self,
                               triangle: CachedPolygon,
                               other: CachedPolygon
                               ) -> Segment[Fraction]:
        """
        Finds the segment shared by two neighbor triangles.

        :param triangle: triangle in the triangular map
        :param other: neighbor of the triangle
        :return: segment that the given triangles share
        """
        return next(edge
                    for edge, neighbor in self.triangular_map[triangle].items()
                    if neighbor == other)

    def unite(self,
              chunk: Chunk,
              other: Chunk) -> Chunk:
        """
        Unites two neighbor chunks into a new one.

        The given chunks are removed from the partition and must not be
        used with it again. The union takes over their neighbors, with
        the edges shared with a neighbor of both added up. Its triangles
        are those of ``chunk`` followed by those of ``other``.
        :param chunk: any chunk from the partition
        :param other: neighbor chunk
        :return: newly created chunk as their union
        :raises ValueError: if the chunks aren't neighbors
        """
        if other not in self.chunk_map[chunk]:
            raise ValueError("Only neighbor chunks can be united")
        union = Chunk(triangles={**chunk.triangles, **other.triangles},
                      area=chunk.area + other.area,
                      perimeter=self.union_perimeter(chunk, other))
        self.chunk_map[union] = self._neighbors_of_union(chunk, other)
        # Outside each other, the united chunks only appear in the
        # entries of the union's neighbors.
        for neighbor, edges in self.chunk_map[union].items():
            neighbors = self.chunk_map[neighbor]
            neighbors.pop(chunk, None)
            neighbors.pop(other, None)
            neighbors[union] = edges
        self.chunk_map.pop(chunk)
        self.chunk_map.pop(other)
        return union

    def union_perimeter(self,
                        chunk: Chunk,
                        other: Chunk) -> float:
        """
        Calculates the perimeter of the union of two distinct chunks.

        The chunks don't have to be neighbors: without shared edges the
        perimeters just add up. Passing the same chunk twice gives a
        wrong result.
        :param chunk: any chunk from the partition
        :param other: any other chunk from the partition
        :return: perimeter of a chunk that could be obtained if the
        chunks were to be united
        """
        shared_length = sum(map(fast_length, self.shared_edges(chunk, other)))
        return chunk.perimeter + other.perimeter - 2 * shared_length

    def _neighbors_of_union(self,
                            chunk: Chunk,
                            other: Chunk) -> SharedEdgesPerNeighbor:
        result = {neighbor: edges
                  for neighbor, edges in self.chunk_map[chunk].items()
                  if neighbor is not other}
        for neighbor, edges in self.chunk_map[other].items():
            if neighbor is not chunk:
                result[neighbor] = result.get(neighbor, frozenset()) | edges
        return result

    def _remove_shared_edge(self,
                            chunk: Chunk,
                            other: Chunk,
                            edge: Segment[Fraction]) -> None:
        self._set_shared_edges(chunk,
                               other,
                               self.shared_edges(chunk, other) - {edge})

    def _add_shared_edge(self,
                         chunk: Chunk,
                         other: Chunk,
                         edge: Segment[Fraction]) -> None:
        self._set_shared_edges(chunk,
                               other,
                               self.shared_edges(chunk, other) | {edge})

    def _set_shared_edges(self,
                          chunk: Chunk,
                          other: Chunk,
                          edges: frozenset[Segment[Fraction]]) -> None:
        # An empty set of shared edges means the chunks stop being
        # neighbors.
        for first, second in ((chunk, other), (other, chunk)):
            if edges:
                self.chunk_map[first][second] = edges
            else:
                self.chunk_map[first].pop(second, None)


    def is_touching_other_chunk(self,
                                triangle: CachedPolygon,
                                *,
                                this_chunk: Chunk) -> bool:
        """
        Checks if the triangle borders a chunk other than its own.

        :param triangle: any triangle from the triangular map
        :param this_chunk: chunk this triangle belongs to
        :return: True if the condition satisfies, False otherwise
        """
        return any(neighbor not in this_chunk.triangles
                   for neighbor in self.triangle_neighbors(triangle))

    def border_triangles(self, chunk: Chunk) -> Iterator[CachedPolygon]:
        """
        Yields the triangles of the chunk that touch other chunks.

        :param chunk: any chunk from the partition
        :return: iterator over triangles on the border of the chunk
        """
        for triangle in chunk.triangles:
            if self.is_touching_other_chunk(triangle, this_chunk=chunk):
                yield triangle

    def movable_triangles(self, chunk: Chunk) -> Iterator[CachedPolygon]:
        """
        Yields the border triangles that can be moved from the chunk.

        :return: iterator over triangles of the given chunk that can be
        moved without causing discontinuity
        """
        yield from (triangle for triangle in self.border_triangles(chunk)
                    if not self.is_articulation_triangle(triangle,
                                                         chunk.triangles))


def to_triangular_map(triangles: list[CachedPolygon]) -> TriangularMap:
    """
    Builds a map of neighboring triangles and their shared segments.

    :param triangles: a list of triangles
    :return: a map that stores relation between neighboring triangles
    and shared segments
    """
    triangular_map: TriangularMap = {}
    triangle_by_edge: dict[Segment[Fraction], CachedPolygon] = {}
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


def to_chunk_map(triangular_map: TriangularMap
                 ) -> dict[Chunk, SharedEdgesPerNeighbor]:
    """
    Converts a triangular map to a chunk map.

    :param triangular_map: input triangular map
    :return: same triangular map but where triangles are wrapped as
    chunks, mapping each neighbor chunk to the set holding the single
    segment shared with it
    """
    chunk_per_triangle = {triangle: to_chunk(triangle)
                          for triangle in triangular_map}
    return {chunk_per_triangle[triangle]:
            {chunk_per_triangle[neighbor]: frozenset([edge])
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
                 perimeter=triangle.perimeter)
