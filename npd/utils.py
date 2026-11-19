"""Various utility functions."""
import math
from collections.abc import (Iterable,
                             Iterator)
from decimal import Decimal
from fractions import Fraction
from functools import reduce
from operator import or_
from typing import TypeVar

from gon.base import (EMPTY,
                      Contour,
                      Empty,
                      Multipolygon,
                      Point,
                      Polygon,
                      Segment)

T = TypeVar("T")


def fast_length(geometry: Segment[Fraction] | Contour[Fraction]) -> float:
    """
    Calculates the length of a segment or a contour using floats.

    gon's ``length`` is exact but symbolic, so slower and, for contours,
    not convertible to float.
    :param geometry: segment or contour to measure
    :return: approximate length
    :raises TypeError: for any other geometry
    """
    if isinstance(geometry, Segment):
        return distance(geometry.start, geometry.end)
    if isinstance(geometry, Contour):
        vertices = geometry.vertices
        return sum(map(distance, vertices, (*vertices[1:], vertices[0])))
    raise TypeError(f"Cannot calculate length of type {type(geometry)}")


def distance(start: Point[Fraction], end: Point[Fraction]) -> float:
    """
    Calculates the distance between two points using floats.

    The differences are exact; only ``math.hypot`` rounds them, which
    keeps it fast for coordinates with many digits.
    :param start: one point
    :param end: another point
    :return: approximate distance
    """
    return math.hypot(start.x - end.x, start.y - end.y)


def inverse_shoelace(area: Fraction,
                     endpoint: Point[Fraction],
                     base_start: Point[Fraction],
                     base_end: Point[Fraction]) -> Point[Fraction]:
    """
    Finds the point on the base forming a triangle of the given area.

    The triangle is formed by the endpoint, the base start and the point.
    Its area grows linearly as the point moves from the base start to
    the base end, up to the area of the triangle formed by the endpoint
    and the whole base, so the point divides the base in the ratio of
    the two areas. It lies strictly inside the base exactly when the
    area is positive and less than that triangle's area; otherwise it
    lies elsewhere on the base's line.
    :param area: area of the triangle to form
    :param endpoint: vertex of the triangle opposite the base
    :param base_start: vertex of the triangle on the base
    :param base_end: other end of the base
    :return: the third vertex of the triangle, on the base's line
    :raises ZeroDivisionError: if the endpoint lies on the base's line
    """
    double_full_area = abs((base_start.x - endpoint.x)
                           * (base_end.y - endpoint.y)
                           - (base_start.y - endpoint.y)
                           * (base_end.x - endpoint.x))
    ratio = 2 * area / double_full_area
    return Point(base_start.x + ratio * (base_end.x - base_start.x),
                 base_start.y + ratio * (base_end.y - base_start.y))


def unite(polygons: Iterable[Polygon[Fraction]]
          ) -> Polygon[Fraction] | Multipolygon[Fraction] | Empty:
    """
    Unites the given polygons into one geometry.

    gon has no n-ary union, so this folds its binary ``|`` operator.

    :param polygons: polygons to unite
    :return: a ``Polygon`` if the union is connected, a ``Multipolygon``
    if it is not, and ``EMPTY`` if no polygons were given
    """
    return reduce(or_, polygons, EMPTY)
