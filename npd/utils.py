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
                      Linear,
                      Multipolygon,
                      Orientation,
                      Point,
                      Polygon,
                      Segment)

T = TypeVar("T")


def fast_length(geometry: Linear[Fraction]) -> float:
    """
    Calculates the length of a segment or a contour using floats.

    :param geometry:
    :return:
    """
    if isinstance(geometry, Segment):
        start = geometry.start
        end = geometry.end
        return math.sqrt((start.x - end.x) ** 2 + (start.y - end.y) ** 2)
    if isinstance(geometry, Contour):
        result = 0.0
        for segment in geometry.segments:
            start = segment.start
            end = segment.end
            result += math.sqrt((start.x - end.x) ** 2
                                + (start.y - end.y) ** 2)
        return result
    raise TypeError(f"Cannot calculate length of type {type(geometry)}")


def to_segments(vertices: list[Point[Fraction]]
                ) -> Iterator[Segment[Fraction]]:
    """Returns an iterator over consecutive segments on the points."""
    yield from map(Segment, vertices, vertices[1:])


def inverse_shoelace(area: Fraction,
                     endpoint: Point[Fraction],
                     base_start: Point[Fraction],
                     base_end: Point[Fraction]) -> Point[Fraction]:
    """
    Finds a point on the base line forming a triangle of the given area.

    The triangle is formed by the endpoint, the base start and the point.
    :param area:
    :param endpoint:
    :param base_start:
    :param base_end:
    :return:
    """
    orientation = Contour([endpoint, base_start, base_end]).orientation
    sign = 1 if orientation is Orientation.COUNTERCLOCKWISE else -1
    if base_end.x - base_start.x == 0:
        numerator = 2 * area + sign * (base_start.x * base_start.y
                                       - endpoint.x * base_start.y)
        denominator = sign * (base_start.x - endpoint.x)
        y_coordinate = numerator / denominator
        return Point(base_start.x, y_coordinate)
    slope, intercept = slope_intercept(base_start, base_end)
    x_diff = base_start.x - endpoint.x
    numerator = 2 * area + sign * (base_start.x * endpoint.y
                                   - endpoint.x * base_start.y
                                   - intercept * x_diff)
    denominator = sign * (endpoint.y - base_start.y + slope * x_diff)
    x_coordinate = numerator / denominator
    y_coordinate = slope * x_coordinate + intercept
    return Point(x_coordinate, y_coordinate)


def slope_intercept(start: Point[Fraction],
                    end: Point[Fraction]) -> tuple[Fraction, Fraction]:
    """
    Calculates the slope and the intercept of the line through two points.

    :param start:
    :param end:
    :return:
    """
    slope = (end.y - start.y) / (end.x - start.x)
    intercept = end.y - slope * end.x
    return slope, intercept


def robust_sqrt(fraction: Fraction) -> Fraction:
    """
    Calculates the square root of a fraction using decimals.

    :param fraction:
    :return:
    """
    return Fraction((fraction.numerator
                     / Decimal(fraction.denominator)).sqrt())


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
