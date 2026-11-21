from fractions import Fraction
from itertools import chain

from gon.base import Polygon
from hypothesis import given

from npd.npd import steiner_points
from tests.strategies.base import (deltas,
                                   varied_polygons)


@given(polygon=varied_polygons,
       delta=deltas)
def test_points_in_polygon(polygon: Polygon[Fraction],
                           delta: Fraction) -> None:
    result = steiner_points(polygon, delta)

    assert all(point in polygon for point in result)


@given(polygon=varied_polygons,
       delta=deltas)
def test_points_distinct(polygon: Polygon[Fraction],
                         delta: Fraction) -> None:
    result = steiner_points(polygon, delta)

    assert len(set(result)) == len(result)


@given(polygon=varied_polygons,
       delta=deltas)
def test_no_polygon_vertices(polygon: Polygon[Fraction],
                             delta: Fraction) -> None:
    vertices = set(chain(polygon.border.vertices,
                         *[hole.vertices for hole in polygon.holes]))

    result = steiner_points(polygon, delta)

    assert vertices.isdisjoint(result)
