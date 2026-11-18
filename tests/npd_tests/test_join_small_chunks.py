from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import (PARTS_COUNT,
                           is_satisfactory_area,
                           join_small_chunks)
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   unit_area_fractions,
                                   varied_polygons)
from tests.utils import (apply,
                         rebuilt_neighbors)

Triangles = frozenset[CachedPolygon]


def to_partition_after(polygon: Polygon[Fraction],
                       delta: Fraction,
                       operations: list[tuple[bool, int, int]]) -> Partition:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    return partition


def triangle_sets(partition: Partition) -> list[Triangles]:
    return [frozenset(chunk.triangles) for chunk in partition.chunk_map]


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       factor=unit_area_fractions)
def test_only_unites_chunks(polygon: Polygon[Fraction],
                            delta: Fraction,
                            operations: list[tuple[bool, int, int]],
                            factor: Fraction) -> None:
    partition = to_partition_after(polygon, delta, operations)
    before = triangle_sets(partition)

    join_small_chunks(partition, requirement=polygon.area * factor)

    after = triangle_sets(partition)
    assert partition.chunk_map == rebuilt_neighbors(partition)
    assert all(sum(triangles <= union for union in after) == 1
               for triangles in before)
    assert len(after) >= min(len(before), PARTS_COUNT)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       factor=unit_area_fractions)
def test_keeps_satisfactory_chunks(polygon: Polygon[Fraction],
                                   delta: Fraction,
                                   operations: list[tuple[bool, int, int]],
                                   factor: Fraction) -> None:
    partition = to_partition_after(polygon, delta, operations)
    requirement = polygon.area * factor
    satisfactory = [frozenset(chunk.triangles)
                    for chunk in partition.chunk_map
                    if is_satisfactory_area(chunk.area,
                                            requirement=requirement)]

    join_small_chunks(partition, requirement=requirement)

    assert all(any(triangles <= chunk.triangles.keys()
                   and is_satisfactory_area(chunk.area,
                                            requirement=requirement)
                   for chunk in partition.chunk_map)
               for triangles in satisfactory)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       factor=unit_area_fractions)
def test_leaves_no_joinable_chunk(polygon: Polygon[Fraction],
                                  delta: Fraction,
                                  operations: list[tuple[bool, int, int]],
                                  factor: Fraction) -> None:
    partition = to_partition_after(polygon, delta, operations)
    requirement = polygon.area * factor

    join_small_chunks(partition, requirement=requirement)

    # Unless it stopped at two chunks, every chunk left unsatisfactory
    # has only neighbors that would overshoot by more than it falls short.
    assert len(partition.chunk_map) <= PARTS_COUNT or all(
        chunk.area + neighbor.area - requirement > requirement - chunk.area
        for chunk in partition.chunk_map
        if not is_satisfactory_area(chunk.area, requirement=requirement)
        for neighbor in partition.chunk_neighbors(chunk))
