from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         is_connected,
                         to_inner_triangle_partition,
                         to_strip_partition,
                         unite_all)


def is_movable(triangle: CachedPolygon,
               chunk: Chunk,
               partition: Partition) -> bool:
    touches_other_chunk = any(
        neighbor is not None and neighbor not in chunk.triangles
        for neighbor in partition.triangular_map[triangle].values())
    rest = [other for other in chunk.triangles if other is not triangle]
    return touches_other_chunk and (not rest
                                    or is_connected(partition, rest))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_movable_triangles(polygon: Polygon[Fraction],
                           delta: Fraction,
                           operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        assert list(partition.movable_triangles(chunk)) == [
            triangle
            for triangle in chunk.triangles
            if is_movable(triangle, chunk, partition)]


def test_single_chunk() -> None:
    partition, chunks = to_inner_triangle_partition()
    chunk = unite_all(partition, *chunks.values())

    assert list(partition.movable_triangles(chunk)) == []


def test_single_triangle() -> None:
    partition, chunks = to_inner_triangle_partition()
    chunk = chunks["inner"]

    assert list(partition.movable_triangles(chunk)) == list(chunk.triangles)


def test_strip() -> None:
    # In each unit square only the triangle touching another square can
    # be moved; the middle square's two triangles both touch one.
    partition, (left, middle, right) = to_strip_partition()

    assert len(list(partition.movable_triangles(left))) == 1
    assert list(partition.movable_triangles(middle)) == list(
        middle.triangles)
    assert len(list(partition.movable_triangles(right))) == 1


def test_articulation_triangle() -> None:
    # `left` touches `inner` but links `top_left` and `top`.
    partition, chunks = to_inner_triangle_partition()
    (left,) = chunks["left"].triangles
    chunk = unite_all(partition,
                      chunks["top_left"], chunks["left"], chunks["top"])

    assert set(partition.movable_triangles(chunk)) == (
        set(chunk.triangles) - {left})
