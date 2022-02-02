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
    context = get_context()
    original_sqrt = context.sqrt
    context._sqrt = math.sqrt
    result = geometry.length
    context._sqrt = original_sqrt
    return result


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


def to_complete_border_vertices(border: Contour,
                                part_vertices: List[Point[Fraction]]
                                ) -> List[Point[Fraction]]:
    """
    Adds points from border of one part of the polygon into the complete
    list of border vertices
    :param border:
    :param part_vertices:
    :return:
    """
    result = list(border.vertices)
    part_vertices_on_border = Multipoint(part_vertices) & border
    extra_vertices = part_vertices_on_border - Multipoint(border.vertices)
    if extra_vertices is EMPTY:
        return list(border.vertices)
    extra_vertices = extra_vertices.points
    if len(extra_vertices) == 2:
        first_extra_vertex, second_extra_vertex = extra_vertices
        first_segment_index, first_segment_start = next(
            ((i, segment.start)
             for i, segment in enumerate(border.segments[1:])
             if first_extra_vertex in segment),
            (len(border.vertices) - 1, border.vertices[-1]))
        second_segment_index = next(
            (i for i, segment in enumerate(border.segments[1:])
             if second_extra_vertex in segment), len(border.vertices) - 1)
        if first_segment_index == second_segment_index:
            first_extra_vertex, second_extra_vertex = sorted(
                (first_extra_vertex, second_extra_vertex),
                key=first_segment_start.distance_to)
            return [*border.vertices[:first_segment_index + 1],
                    first_extra_vertex, second_extra_vertex,
                    *border.vertices[first_segment_index + 1:]]
        if first_segment_index < second_segment_index:
            result.insert(second_segment_index + 1, second_extra_vertex)
            result.insert(first_segment_index + 1, first_extra_vertex)
        else:
            result.insert(first_segment_index + 1, first_extra_vertex)
            result.insert(second_segment_index + 1, second_extra_vertex)
        return result
    extra_vertex = extra_vertices.pop()
    i = next((i for i, segment in enumerate(border.segments[1:])
              if extra_vertex in segment), len(border.vertices) - 1)
    return [*border.vertices[:i + 1], extra_vertex, *border.vertices[i + 1:]]


def get_splitter_vertices(contour: Contour[Fraction],
                          other: Contour[Fraction]) -> List[Point[Fraction]]:
    """
    :param contour: contour delimiting one part of the polygon
    :param other: contour delimiting the other part of the polygon
    :return: ordered list of vertices that are shared between the
    two parts of the polygon
    """
    shared_border = contour & other
    shared_segments = (list(shared_border.segments)
                       if isinstance(shared_border, Multisegment)
                       else [shared_border])
    unordered_splitter_vertices = set(chain.from_iterable(
        (segment.start, segment.end) for segment in shared_segments))
    start = next((i for i, v in enumerate(contour.vertices)
                  if v not in unordered_splitter_vertices), None)
    if start is None:
        own_part = contour - shared_border
        j = next(i for i, v in enumerate(contour.vertices) if v in own_part)
        start = (j if contour.vertices[j - 1] in own_part
                 else (j + 1) % len(contour.vertices))
    return list(arrange_along(unordered_splitter_vertices,
                              contour,
                              start=start))


def arrange_along(vertices: Set[Point[Fraction]],
                  contour: Contour[Fraction],
                  *,
                  start: int) -> Iterator[Point[Fraction]]:
    """
    Arranges given vertices in the order as they are encountered along
    the contour in the counterclockwise order starting from the vertex
    with the given index
    :param vertices: vertices to be arranged
    :param contour: input contour
    :param start: index of a vertex on a counter from where the
    arrangement of the vertices starts
    :return: arranged vertices
    """
    contour_vertices = shift(list(contour.vertices), start)
    contour_vertices.append(contour_vertices[0])
    segments = to_segments(contour_vertices)
    for segment in segments:
        vertices_on_segment = {vertex for vertex in vertices
                               if vertex in segment}
        vertices -= vertices_on_segment
        ordered_vertices_on_segment = sorted(vertices_on_segment,
                                             key=segment.start.distance_to)
        yield from iter(ordered_vertices_on_segment)


def compactness(*, area_: float, perimeter: float) -> float:
    """
    :param area_:
    :param perimeter:
    :return:
    """
    return area_ ** 0.5 / perimeter


def unite(geometries: Iterable[Polygon[Fraction]]) -> Polygon[Fraction]:
    """
    :param geometries:
    :return:
    """
    return reduce(or_, geometries, EMPTY)


def shift(values: List[T], step: int) -> List[T]:
    """
    Shifts the values in the given list by the given step
    :param values: input values
    :param step: indicates how many positions will be shifted
    :return: a list with shifted values
    """
    return values[step:] + values[:step] if step else values
