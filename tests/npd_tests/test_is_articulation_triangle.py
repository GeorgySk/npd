from fractions import Fraction

import pytest
from gon.base import (Contour,
                      Point,
                      Polygon)
from hypothesis import given

from npd.npd import to_partition
from npd.structures import Partition
from tests.strategies.base import (deltas,
                                   partition_operations,
                                   polygons,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         is_connected)


@given(polygon=polygons,
       delta=deltas)
def test_non_articulation_node_presence(polygon: Polygon[Fraction],
                                        delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)
    # There is no need to split partition into chunks and using only a
    # particular chunk for the test.
    assert any(
        not partition.is_articulation_triangle(
            triangle, chunk_triangles=partition.triangular_map)
        for triangle in partition.triangular_map)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_matches_connectivity(polygon: Polygon[Fraction],
                              delta: Fraction,
                              operations: list[tuple[bool, int, int]]
                              ) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        if len(chunk.triangles) < 2:  # noqa: PLR2004
            continue
        for triangle in chunk.triangles:
            rest = {other: None
                    for other in chunk.triangles
                    if other is not triangle}
            assert (partition.is_articulation_triangle(triangle,
                                                       chunk.triangles)
                    is not is_connected(partition, rest))


def to_single_triangle_partition() -> Partition:
    border = Contour([Point(Fraction(0), Fraction(0)),
                      Point(Fraction(1), Fraction(0)),
                      Point(Fraction(0), Fraction(1))])
    return Partition(Polygon(border, []))


def test_single_triangle() -> None:
    partition = to_single_triangle_partition()
    (triangle,) = partition.triangular_map

    assert not partition.is_articulation_triangle(triangle,
                                                  {triangle: None})


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_rejects_foreign_triangle(polygon: Polygon[Fraction],
                                  delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)
    triangles = list(partition.triangular_map)

    with pytest.raises(ValueError, match="one of the given"):
        partition.is_articulation_triangle(triangles[0], triangles[1:])
