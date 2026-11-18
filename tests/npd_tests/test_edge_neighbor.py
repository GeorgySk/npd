from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import edge_neighbor
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


def touches_other_chunk(triangle: CachedPolygon,
                        chunk: Chunk,
                        partition: Partition) -> bool:
    return any(other is not None and other not in chunk.triangles
               for other in partition.triangular_map[triangle].values())


def meets_condition(candidate: CachedPolygon,
                    triangle: CachedPolygon,
                    chunk: Chunk,
                    partition: Partition) -> bool:
    """The condition as stated in the docstring."""
    return not any(
        touches_other_chunk(other, chunk, partition)
        for other in partition.triangular_map[candidate].values()
        if other is not None and other is not triangle)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_returns_first_qualifying_neighbor(polygon: Polygon[Fraction],
                                           delta: Fraction,
                                           operations: list[tuple[bool,
                                                                  int,
                                                                  int]]
                                           ) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        for triangle in chunk.triangles:
            if not touches_other_chunk(triangle, chunk, partition):
                continue
            expected = next(
                (candidate
                 for candidate in partition.triangle_neighbors(triangle)
                 if candidate in chunk.triangles
                 and meets_condition(candidate, triangle, chunk, partition)),
                None)
            if expected is None:
                with pytest.raises(ValueError, match="Couldn't find"):
                    edge_neighbor(triangle, chunk, partition)
            else:
                assert edge_neighbor(triangle, chunk, partition) is expected


def test_skips_neighbor_near_border() -> None:
    partition, chunks = to_inner_triangle_partition()
    (bottom,) = chunks["bottom"].triangles
    (bottom_left,) = chunks["bottom_left"].triangles
    chunk = unite_all(partition, *(chunk for name, chunk in chunks.items()
                                   if name != "bottom_right"))

    # `inner` is skipped: its neighbor `right` touches `bottom_right`.
    assert edge_neighbor(bottom, chunk, partition) is bottom_left


def test_no_qualifying_neighbor() -> None:
    partition, chunks = to_inner_triangle_partition()
    (inner,) = chunks["inner"].triangles
    chunk = unite_all(partition,
                      chunks["inner"], chunks["bottom"], chunks["left"])

    with pytest.raises(ValueError, match="Couldn't find"):
        edge_neighbor(inner, chunk, partition)
