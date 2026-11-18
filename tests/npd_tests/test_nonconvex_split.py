from fractions import Fraction

import pytest
from gon.base import (Contour,
                      Point,
                      Polygon)
from hypothesis import given

from npd import nonconvex
from npd.npd import to_partition
from tests.utils import (is_connected,
                         rebuilt_neighbors)
from tests.strategies.base import (high_requirement_factors,
                                   invalid_requirement_factors,
                                   polygons,
                                   sparse_deltas,
                                   varied_polygons)


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


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       factor=high_requirement_factors)
def test_high_requirement(polygon: Polygon[Fraction],
                          delta: Fraction,
                          factor: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)
    if len(partition.chunk_map) < 2:  # noqa: PLR2004
        return  # a single triangle can't be split yet
    requirement = polygon.area * factor

    chunk, remainder = nonconvex.split(partition, requirement=requirement)

    assert chunk.area == requirement
    assert remainder.area == polygon.area - requirement


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
