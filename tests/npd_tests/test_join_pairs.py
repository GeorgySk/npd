from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import join_pairs
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import apply

Triangles = frozenset[CachedPolygon]


def to_partition_after(polygon: Polygon[Fraction],
                       delta: Fraction,
                       operations: list[tuple[bool, int, int]]) -> Partition:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    return partition


def snapshot(partition: Partition) -> dict[Chunk, Triangles]:
    return {chunk: frozenset(chunk.triangles) for chunk in partition.chunk_map}


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_merges_pairs_of_neighbors(polygon: Polygon[Fraction],
                                   delta: Fraction,
                                   operations: list[tuple[bool, int, int]]
                                   ) -> None:
    partition = to_partition_after(polygon, delta, operations)
    before = snapshot(partition)
    smallest_size = min(map(len, before.values()))
    pairs = {before[chunk] | before[neighbor]: (chunk, neighbor)
             for chunk in before
             for neighbor in partition.chunk_neighbors(chunk)}

    join_pairs(partition)

    for chunk, triangles in snapshot(partition).items():
        if chunk in before:
            assert triangles == before[chunk]
        else:
            first, second = pairs[triangles]
            assert smallest_size in (len(before[first]), len(before[second]))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_unmerged_sources_have_no_unmerged_neighbors(
        polygon: Polygon[Fraction],
        delta: Fraction,
        operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition_after(polygon, delta, operations)
    before = snapshot(partition)
    smallest_size = min(map(len, before.values()))
    neighbors = {chunk: partition.chunk_neighbors(chunk) for chunk in before}

    join_pairs(partition)

    assert all(neighbor not in partition.chunk_map
               for chunk in before
               if chunk in partition.chunk_map
               and len(before[chunk]) == smallest_size
               for neighbor in neighbors[chunk])
