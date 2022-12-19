"""Various utility functions"""
import math
from decimal import Decimal
from fractions import Fraction
from functools import reduce
from itertools import chain
from operator import or_
from typing import (Iterable,
                    Iterator,
                    List,
                    Set,
                    Tuple,
                    TypeVar)

from gon.base import (Contour,
                      EMPTY,
                      Linear,
                      Multipoint,
                      Multisegment,
                      Orientation,
                      Point,
                      Polygon,
                      Segment)
from ground.base import get_context

T = TypeVar('T')


def fast_length(geometry: Linear[Fraction]) -> float:
    """
    :param geometry:
    :return:
    """
    if isinstance(geometry, Segment):
        start = geometry.start
        end = geometry.end
        return math.sqrt((start.x - end.x) ** 2 + (start.y - end.y) ** 2)
    elif isinstance(geometry, Contour):
        result = 0
        for segment in geometry.segments:
            start = segment.start
            end = segment.end
            result += math.sqrt((start.x - end.x) ** 2
                                + (start.y - end.y) ** 2)
        return result
    else:
        raise TypeError(f"Cannot calculate length of type {type(geometry)}")
    # context = get_context()
    # original_sqrt = context.sqrt
    # context._sqrt = math.sqrt
    # result = geometry.length
    # context._sqrt = original_sqrt
    # return result


def to_segments(vertices: List[Point]) -> Iterator[Segment]:
    """
    Returns an iterator over consecutive segments
    built on the given list of points
    """
    yield from map(Segment, vertices, vertices[1:])


def inverse_shoelace(area: Fraction,
                     endpoint: Point,
                     base_start: Point,
                     base_end: Point) -> Point:
    """
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


def slope_intercept(start: Point, end: Point) -> Tuple[Fraction, Fraction]:
    """
    :param start:
    :param end:
    :return:
    """
    slope = (end.y - start.y) / (end.x - start.x)
    intercept = end.y - slope * end.x
    return slope, intercept


def robust_sqrt(fraction: Fraction) -> Fraction:
    """
    :param fraction:
    :return:
    """
    return Fraction((fraction.numerator
                     / Decimal(fraction.denominator)).sqrt())


def unite(geometries: Iterable[Polygon[Fraction]]) -> Polygon[Fraction]:
    """
    :param geometries:
    :return:
    """
    return reduce(or_, geometries, EMPTY)
