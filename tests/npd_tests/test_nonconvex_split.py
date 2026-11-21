import math
from fractions import Fraction

import pytest
from gon.base import (Contour,
                      Point,
                      Polygon)
from hypothesis import given

from npd import nonconvex
from npd.npd import to_partition
from npd.structures import (Chunk,
                            Partition,
                            to_triangular_map)
from tests.strategies.base import (high_requirement_factors,
                                   invalid_requirement_factors,
                                   polygons,
                                   requirement_factors,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (border_length,
                         is_connected,
                         rebuilt_neighbors)


@given(polygon=polygons,
       delta=sparse_deltas,
       factor=invalid_requirement_factors)
def test_invalid_requirement(polygon: Polygon[Fraction],
                             delta: Fraction,
                             factor: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    with pytest.raises(ValueError, match="positive and less than"):
        nonconvex.split(partition, requirement=polygon.area * factor)


def test_requirement_leaving_two_unsatisfactory_chunks() -> None:
    # Joining reaches two chunks, neither of them satisfactory; uniting
    # them would leave a single chunk to split.
    square = Polygon(Contour([Point(Fraction(0), Fraction(0)),
                              Point(Fraction(2), Fraction(0)),
                              Point(Fraction(2), Fraction(2)),
                              Point(Fraction(0), Fraction(2))]), [])
    partition = to_partition(square, delta=Fraction(1))
    requirement = square.area * Fraction(19, 20)

    chunk, remainder = nonconvex.split(partition, requirement=requirement)

    assert chunk.area == requirement
    assert remainder.area == square.area - requirement


def assert_valid_split(partition: Partition,
                       chunk: Chunk,
                       remainder: Chunk,
                       requirement: Fraction) -> None:
    assert chunk.area == requirement
    assert remainder.area == partition.total_area - requirement
    assert set(partition.chunk_map) == {chunk, remainder}
    assert is_connected(partition, chunk.triangles)
    assert is_connected(partition, remainder.triangles)
    assert partition.chunk_map == rebuilt_neighbors(partition)
    assert partition.triangular_map == to_triangular_map(
        list(partition.triangular_map))
    for chunk_ in (chunk, remainder):
        assert chunk_.area == sum(triangle.area
                                  for triangle in chunk_.triangles)
        assert math.isclose(chunk_.perimeter,
                            border_length(chunk_.triangles, partition))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       factor=requirement_factors)
def test_any_requirement(polygon: Polygon[Fraction],
                         delta: Fraction,
                         factor: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)
    if len(partition.triangular_map) < 2:  # noqa: PLR2004
        return
    requirement = polygon.area * factor

    chunk, remainder = nonconvex.split(partition, requirement=requirement)

    assert_valid_split(partition, chunk, remainder, requirement)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       factor=high_requirement_factors)
def test_high_requirement(polygon: Polygon[Fraction],
                          delta: Fraction,
                          factor: Fraction) -> None:
    # Joining can end at two chunks, neither of them satisfactory.
    partition = to_partition(polygon, delta=delta)
    if len(partition.triangular_map) < 2:  # noqa: PLR2004
        return
    requirement = polygon.area * factor

    chunk, remainder = nonconvex.split(partition, requirement=requirement)

    assert_valid_split(partition, chunk, remainder, requirement)


def test_single_triangle() -> None:
    triangle = Polygon(Contour([Point(Fraction(0), Fraction(0)),
                                Point(Fraction(1), Fraction(0)),
                                Point(Fraction(0), Fraction(1))]), [])
    partition = Partition(triangle)

    with pytest.raises(ValueError, match="at least two triangles"):
        nonconvex.split(partition, requirement=triangle.area / 2)


def test_l_shape() -> None:
    # `readjust` moves triangles it had set aside to try last; they used
    # to be checked as articulation triangles of a chunk they had left.
    l_shape = Polygon(Contour([Point(Fraction(x), Fraction(y))
                               for x, y in [(0, 0), (2, 0), (2, 1),
                                            (1, 1), (1, 2), (0, 2)]]), [])
    partition = to_partition(l_shape, delta=Fraction(1, 3))
    requirement = Fraction(1, 2)

    chunk, remainder = nonconvex.split(partition, requirement=requirement)

    assert chunk.area == requirement
    assert remainder.area == l_shape.area - requirement
    assert set(partition.chunk_map) == {chunk, remainder}
    assert set(chunk.triangles) | set(remainder.triangles) == set(
        partition.triangular_map)
    assert is_connected(partition, chunk.triangles)
    assert is_connected(partition, remainder.triangles)


def test_l_shape_keeps_chunk_map() -> None:
    # The final cut also cuts a triangle of the remainder, splitting an
    # edge the chunks share.
    l_shape = Polygon(Contour([Point(Fraction(x), Fraction(y))
                               for x, y in [(0, 0), (2, 0), (2, 1),
                                            (1, 1), (1, 2), (0, 2)]]), [])
    partition = to_partition(l_shape, delta=Fraction(1, 4))

    nonconvex.split(partition, requirement=Fraction(7, 10))

    assert partition.chunk_map == rebuilt_neighbors(partition)


def test_pentagon_needing_triangle_split() -> None:
    # `add_triangles` runs out of movable triangles and has to split one,
    # but `edge_neighbor` used to find no neighbor to split it with.
    coordinates: list[tuple[int | Fraction, int]] = [
        (0, 0), (2, 1), (3, 4), (1, 4), (Fraction(3, 2), 1)
    ]
    pentagon = Polygon(Contour([Point(Fraction(x), Fraction(y))
                                for x, y in coordinates]),
                       [])
    partition = to_partition(pentagon, delta=Fraction(52, 59))
    requirement = Fraction(7, 2)

    chunk, remainder = nonconvex.split(partition, requirement=requirement)

    assert chunk.area == requirement
    assert remainder.area == pentagon.area - requirement
    assert is_connected(partition, chunk.triangles)
    assert is_connected(partition, remainder.triangles)
