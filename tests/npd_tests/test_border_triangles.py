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
                         to_inner_triangle_partition,
                         unite_all)


def to_owners(partition: Partition) -> dict[CachedPolygon, Chunk]:
    return {triangle: chunk
            for chunk in partition.chunk_map
            for triangle in chunk.triangles}


def touches_other_chunk(triangle: CachedPolygon,
                        owners: dict[CachedPolygon, Chunk],
                        partition: Partition) -> bool:
    return any(neighbor is not None
               and owners[neighbor] is not owners[triangle]
               for neighbor in partition.triangular_map[triangle].values())


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_is_touching_other_chunk(polygon: Polygon[Fraction],
                                 delta: Fraction,
                                 operations: list[tuple[bool, int, int]]
                                 ) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    owners = to_owners(partition)

    for triangle, chunk in owners.items():
        assert partition.is_touching_other_chunk(
            triangle, this_chunk=chunk) is touches_other_chunk(
            triangle, owners, partition)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_border_triangles(polygon: Polygon[Fraction],
                          delta: Fraction,
                          operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    owners = to_owners(partition)

    for chunk in partition.chunk_map:
        assert list(partition.border_triangles(chunk)) == [
            triangle
            for triangle in chunk.triangles
            if touches_other_chunk(triangle, owners, partition)]


def test_single_chunk() -> None:
    partition, chunks = to_inner_triangle_partition()
    chunk = unite_all(partition, *chunks.values())

    assert list(partition.border_triangles(chunk)) == []
    assert not any(partition.is_touching_other_chunk(triangle,
                                                     this_chunk=chunk)
                   for triangle in chunk.triangles)


def test_inner_triangles() -> None:
    # Uniting all but `inner` leaves `bottom`, `right` and `left` as the
    # only triangles touching it.
    partition, chunks = to_inner_triangle_partition()
    inner = chunks.pop("inner")
    chunk = unite_all(partition, *chunks.values())

    assert set(partition.border_triangles(chunk)) == {
        triangle
        for name in ("bottom", "right", "left")
        for triangle in chunks[name].triangles}
    assert list(partition.border_triangles(inner)) == list(inner.triangles)
