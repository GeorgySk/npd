import math
from collections.abc import Collection
from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon,
                      Segment)

from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)


def apply(partition: Partition,
          operations: list[tuple[bool, int, int]]) -> None:
    for is_move, source_index, choice_index in operations:
        if is_move:
            owners = {triangle: chunk
                      for chunk in partition.chunk_map
                      for triangle in chunk.triangles}
            # a move must leave its source non-empty
            moves = [(triangle, source, owners[neighbor])
                     for source in partition.chunk_map
                     if len(source.triangles) > 1
                     for triangle in partition.movable_triangles(source)
                     for neighbor in partition.triangle_neighbors(triangle)
                     if owners[neighbor] is not source]
            if not moves:
                continue
            triangle, source, target = moves[choice_index % len(moves)]
            partition.move(triangle, source=source, target=target)
        else:
            chunks = list(partition.chunk_map)
            source = chunks[source_index % len(chunks)]
            neighbors = partition.chunk_neighbors(source)
            if not neighbors:
                continue
            partition.unite(source,
                            neighbors[choice_index % len(neighbors)])


def rebuilt_neighbors(partition: Partition
                      ) -> dict[Chunk, dict[Chunk,
                                            frozenset[Segment[Fraction]]]]:
    owners: dict[CachedPolygon, Chunk] = {
        triangle: chunk
        for chunk in partition.chunk_map
        for triangle in chunk.triangles}
    result: dict[Chunk, dict[Chunk, set[Segment[Fraction]]]] = {
        chunk: {} for chunk in partition.chunk_map}
    for triangle, neighbors in partition.triangular_map.items():
        chunk = owners[triangle]
        for edge, neighbor in neighbors.items():
            if neighbor is not None and owners[neighbor] is not chunk:
                result[chunk].setdefault(owners[neighbor], set()).add(edge)
    return {chunk: {neighbor: frozenset(edges)
                    for neighbor, edges in neighbors_edges.items()}
            for chunk, neighbors_edges in result.items()}


def is_connected(partition: Partition,
                 triangles: Collection[CachedPolygon]) -> bool:
    """Checks if the triangles are connected through shared edges."""
    root = next(iter(triangles))
    seen, stack = {root}, [root]
    while stack:
        for neighbor in partition.triangle_neighbors(stack.pop()):
            if neighbor in triangles and neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == len(triangles)


def border_length(triangles: Collection[CachedPolygon],
                  partition: Partition) -> float:
    return sum(math.dist((edge.start.x, edge.start.y),
                         (edge.end.x, edge.end.y))
               for triangle in triangles
               for edge, neighbor in partition.triangular_map[triangle].items()
               if neighbor not in triangles)


def to_strip_partition() -> tuple[Partition, list[Chunk]]:
    """
    Builds a 3 x 1 strip cut into unit squares, one chunk per square.

    Only the middle chunk is an articulation chunk.
    """
    border = Contour([Point(Fraction(0), Fraction(0)),
                      Point(Fraction(3), Fraction(0)),
                      Point(Fraction(3), Fraction(1)),
                      Point(Fraction(0), Fraction(1))])
    partition = Partition(Polygon(border, []),
                          extra_points=[Point(Fraction(x), Fraction(y))
                                        for x in (1, 2)
                                        for y in (0, 1)])
    squares: dict[Fraction, list[Chunk]] = {}
    for chunk in partition.chunk_map:
        (triangle,) = chunk.triangles
        left_x = min(vertex.x for vertex in triangle.border.vertices)
        squares.setdefault(left_x, []).append(chunk)
    return partition, [partition.unite(*squares[Fraction(x)])
                       for x in range(3)]


def to_point(x: int, y: int) -> Point[Fraction]:
    return Point(Fraction(x), Fraction(y))


def to_inner_triangle_partition() -> tuple[Partition,
                                           dict[str, Chunk]]:
    """
    Square with an inner triangle, one chunk per triangle.

    The inner triangle ``inner`` borders ``bottom``, ``right`` and
    ``left``, which have no edge on the square's border; each of the
    four corner triangles borders two of them or another corner.
    """
    partition = Partition(
        Polygon(Contour([to_point(0, 0), to_point(4, 0),
                         to_point(4, 4), to_point(0, 4)]), []),
        extra_points=[to_point(1, 1), to_point(3, 1), to_point(2, 3)])
    names = {
        frozenset({(1, 1), (3, 1), (2, 3)}): "inner",
        frozenset({(1, 1), (4, 0), (3, 1)}): "bottom",
        frozenset({(2, 3), (3, 1), (4, 4)}): "right",
        frozenset({(0, 4), (1, 1), (2, 3)}): "left",
        frozenset({(3, 1), (4, 0), (4, 4)}): "bottom_right",
        frozenset({(0, 4), (2, 3), (4, 4)}): "top",
        frozenset({(0, 0), (1, 1), (0, 4)}): "top_left",
        frozenset({(0, 0), (4, 0), (1, 1)}): "bottom_left",
    }
    chunks = {}
    for chunk in partition.chunk_map:
        (triangle,) = chunk.triangles
        vertices = frozenset((int(vertex.x), int(vertex.y))
                             for vertex in triangle.border.vertices)
        chunks[names[vertices]] = chunk
    return partition, chunks


def unite_all(partition: Partition, *chunks: Chunk) -> Chunk:
    """Unites connected chunks, each time with a neighbor of the union."""
    result, *rest = chunks
    while rest:
        chunk = next(chunk for chunk in rest
                     if chunk in partition.chunk_map[result])
        rest.remove(chunk)
        result = partition.unite(result, chunk)
    return result
