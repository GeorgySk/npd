import math
from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon)
from hypothesis import given

from npd.nonconvex import swap_leaves
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         rebuilt_neighbors)

# moved triangle, source, target, and the triangle's neighbors in total,
# in the source and in the target at the time of the move
Move = tuple[CachedPolygon, Chunk, Chunk, int, int, int]


def to_partition_after(polygon: Polygon[Fraction],
                       delta: Fraction,
                       operations: list[tuple[bool, int, int]]) -> Partition:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    return partition


def swap_recording_moves(partition: Partition) -> list[Move]:
    moves: list[Move] = []
    move = partition.move

    def recording_move(triangle: CachedPolygon,
                       *,
                       source: Chunk,
                       target: Chunk) -> None:
        neighbors = partition.triangle_neighbors(triangle)
        moves.append((triangle,
                      source,
                      target,
                      len(neighbors),
                      sum(neighbor in source.triangles
                          for neighbor in neighbors),
                      sum(neighbor in target.triangles
                          for neighbor in neighbors)))
        move(triangle, source=source, target=target)

    partition.move = recording_move  # type: ignore[method-assign]
    swap_leaves(partition)
    return moves


def is_connected(partition: Partition, chunk: Chunk) -> bool:
    root = next(iter(chunk.triangles))
    seen, stack = {root}, [root]
    while stack:
        for neighbor in partition.triangle_neighbors(stack.pop()):
            if neighbor in chunk.triangles and neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == len(chunk.triangles)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_moves_leaves(polygon: Polygon[Fraction],
                      delta: Fraction,
                      operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition_after(polygon, delta, operations)

    moves = swap_recording_moves(partition)

    assert all((total, in_source, in_target) == (3, 1, 2)
               for _, _, _, total, in_source, in_target in moves)
    sources = [source for _, source, *_ in moves]
    assert len(sources) == len(set(sources))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_keeps_chunks_valid(polygon: Polygon[Fraction],
                            delta: Fraction,
                            operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition_after(polygon, delta, operations)
    connected = [chunk for chunk in partition.chunk_map
                 if is_connected(partition, chunk)]

    swap_leaves(partition)

    assert partition.chunk_map == rebuilt_neighbors(partition)
    assert all(chunk.triangles for chunk in partition.chunk_map)
    assert all(is_connected(partition, chunk) for chunk in connected)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_perimeters_dont_grow(polygon: Polygon[Fraction],
                              delta: Fraction,
                              operations: list[tuple[bool, int, int]]
                              ) -> None:
    partition = to_partition_after(polygon, delta, operations)
    before = {chunk: chunk.perimeter for chunk in partition.chunk_map}

    swap_leaves(partition)

    assert all(chunk.perimeter <= before[chunk]
               or math.isclose(chunk.perimeter, before[chunk])
               for chunk in partition.chunk_map)


def test_swaps_inner_leaf() -> None:
    def to_point(x: int, y: int) -> Point[Fraction]:
        return Point(Fraction(x), Fraction(y))

    # The extra points form an inner triangle, surrounded by three
    # triangles with no edge on the border.
    partition = Partition(
        Polygon(Contour([to_point(0, 0), to_point(4, 0),
                         to_point(4, 4), to_point(0, 4)]), []),
        extra_points=[to_point(1, 1), to_point(3, 1), to_point(2, 3)])
    inner = next(triangle for triangle in partition.triangular_map
                 if set(triangle.border.vertices)
                 == {to_point(1, 1), to_point(3, 1), to_point(2, 3)})
    *outer_neighbors, other_neighbor = partition.triangle_neighbors(inner)
    owners = {next(iter(chunk.triangles)): chunk
              for chunk in partition.chunk_map}
    # The inner triangle and one neighbor in one chunk, the rest in another.
    leaf_chunk = partition.unite(owners[inner], owners[other_neighbor])
    rest = owners[outer_neighbors[0]]
    while set(partition.chunk_map) != {leaf_chunk, rest}:
        rest = partition.unite(rest, next(
            chunk for chunk in partition.chunk_neighbors(rest)
            if chunk is not leaf_chunk))
    assert {inner, other_neighbor} == set(leaf_chunk.triangles)

    swap_leaves(partition)

    # Both triangles of the leaf chunk are leaves; one of them moves.
    assert len(leaf_chunk.triangles) == 1
    assert set(rest.triangles) == (set(partition.triangular_map)
                                   - set(leaf_chunk.triangles))
