from fractions import Fraction

import pytest
from gon.base import (Point,
                      Polygon)
from hypothesis import given

from npd.utils import touching_points
from tests.strategies.base import polygons
from tests.utils import (HOLED_SQUARE,
                         TOUCHING_HOLES_POLYGONS,
                         to_point,
                         to_polygon)


@given(polygons)
def test_no_holes(polygon: Polygon[Fraction]) -> None:
    assert touching_points(polygon) == []


@pytest.mark.parametrize(("polygon", "point"), TOUCHING_HOLES_POLYGONS)
def test_touching_inside_edge(polygon: Polygon[Fraction],
                              point: Point[Fraction]) -> None:
    assert touching_points(polygon) == [point]


def test_touching_several_contours() -> None:
    # Both holes touch the border's bottom edge, and each other, at (2, 0).
    polygon = to_polygon([(0, 0), (4, 0), (4, 4), (0, 4)],
                         [(2, 0), (1, 1), (1, 2)],
                         [(2, 0), (3, 2), (3, 1)])

    assert touching_points(polygon) == [to_point(2, 0)]


@pytest.mark.parametrize("polygon", [
    HOLED_SQUARE,
    # a hole touching the border at a vertex of both
    to_polygon([(0, 0), (3, 0), (3, 3), (0, 3)],
               [(1, 1), (1, 2), (3, 3)]),
    # holes touching each other at a vertex of both
    to_polygon([(0, 0), (4, 0), (4, 4), (0, 4)],
               [(1, 1), (1, 2), (2, Fraction(3, 2))],
               [(2, Fraction(3, 2)), (3, 1), (3, 2)])])
def test_not_touching_inside_edge(polygon: Polygon[Fraction]) -> None:
    assert touching_points(polygon) == []
