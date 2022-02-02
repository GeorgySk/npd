"""Smoothing of the border between two polygons"""
import logging
from fractions import Fraction
from functools import partial
from typing import (List,
                    Tuple)

from gon.base import (Contour,
                      Orientation,
                      Point,
                      Polygon,
                      Segment)
from ground.base import get_context

from draw import draw
from npd import convex

LOGGER = logging.getLogger(__name__)
context = get_context()


def smooth(splitter_vertices: List[Point[Fraction]],
           part: Polygon,
           other: Polygon,
           border_vertices: List[Point[Fraction]],
           polygon: Polygon) -> Tuple[Polygon, Polygon]:
    """
    :param splitter_vertices:
    :param part:
    :param other:
    :param border_vertices:
    :param polygon:
    :return:
    """
    splitter_vertices_on_holes = {vertex for vertex in splitter_vertices
                                  if any(vertex in hole
                                         for hole in polygon.holes)}
    original_part = part
    LOGGER.debug("Start of a smoothing process.")
    idx = 1  # index of a point to remove
    draw(part, other, red=splitter_vertices[0])
    while True:
        LOGGER.debug("New iteration. Index: %s", idx)
        draw(part, other, splitter_vertices)
        draw(part, other, border_vertices)
        if splitter_vertices[idx] in border_vertices:
            if part == original_part or all(v in border_vertices
                                            for v in splitter_vertices):
                LOGGER.debug("Smoothing process finished.")
                return part, other
            LOGGER.debug("Repeat smoothing procedure.")
            original_part = part
            idx = 2
        if splitter_vertices[idx] in splitter_vertices_on_holes:
            LOGGER.debug("Vertex is on a hole. Skipping...")
            idx += 1
            continue
        triangle_contour = Contour(splitter_vertices[idx - 1:idx + 2])
        draw(part, other, red=triangle_contour)
        orientation = triangle_contour.orientation
        area = Polygon(triangle_contour, []).area
        LOGGER.debug("Found triangle to be removed with %s and area %s.",
                     str(orientation),
                     area)
        if splitter_vertices[idx - 1] in border_vertices:
            if orientation is Orientation.COUNTERCLOCKWISE:
                to_prepend = border_vertices[(border_vertices.index(
                    (splitter_vertices[idx - 1])) + 1) % len(border_vertices)]
            else:
                to_prepend = border_vertices[(border_vertices.index(
                    splitter_vertices[idx - 1])) - 1]
            if idx == 1:
                splitter_vertices.insert(0, to_prepend)
                idx += 1
            elif idx == 2:
                splitter_vertices[0] = to_prepend
            else:
                raise RuntimeError
        if splitter_vertices[idx + 1] in border_vertices:
            if orientation is Orientation.COUNTERCLOCKWISE:
                to_append = border_vertices[(border_vertices.index(
                    splitter_vertices[idx + 1])) - 1]
            else:
                to_append = border_vertices[((border_vertices.index(
                    splitter_vertices[idx + 1])) + 1) % len(border_vertices)]
            if len(splitter_vertices) == idx + 3:
                splitter_vertices[-1] = to_append
            elif len(splitter_vertices) == idx + 2:
                splitter_vertices.append(to_append)
            else:
                raise RuntimeError
        tail_segment, head_segment, area_diff = (
            to_ccw_segments(area=area,
                            border_vertices=border_vertices,
                            tail_start=splitter_vertices[idx - 1],
                            tail_vertex=splitter_vertices[idx - 2],
                            head_end=splitter_vertices[idx + 1],
                            head_vertex=splitter_vertices[idx + 2])
            if orientation is Orientation.COUNTERCLOCKWISE
            else to_cw_segments(area=area,
                                tail_end=splitter_vertices[idx - 1],
                                tail_vertex=splitter_vertices[idx - 2],
                                head_start=splitter_vertices[idx + 1],
                                head_vertex=splitter_vertices[idx + 2]))
        if any(vertex not in polygon
               for vertex in [head_segment.start, head_segment.end,
                              tail_segment.start, tail_segment.end]):
            LOGGER.debug("Segments cross polygon border. Skipping...")
            # Is there a better option than to skip?
            idx += 1
            continue
        draw(part, other, red=triangle_contour,
             blue=[head_segment, tail_segment])
        domain, countersegment = (
            (head_segment, tail_segment)
            if orientation is Orientation.COUNTERCLOCKWISE
            else (tail_segment, head_segment))
        splitter = convex.to_splitter(domain, countersegment, area - area_diff)
        if orientation is Orientation.CLOCKWISE:
            splitter = Segment(splitter.end, splitter.start)
        LOGGER.debug("Splitter line is found.")
        draw(part, other, blue=[head_segment, tail_segment], red=splitter)
        part = update_part(part,
                           splitter_vertices=splitter_vertices,
                           idx=idx,
                           polygon_border_vertices=polygon.border.vertices,
                           border_vertices=border_vertices,
                           splitter=splitter)
        part.validate()
        assert part.area == original_part.area
        LOGGER.debug("First part of the polygon is updated.")
        other = update_part(other,
                            splitter_vertices=splitter_vertices,
                            idx=idx,
                            polygon_border_vertices=polygon.border.vertices,
                            border_vertices=border_vertices,
                            splitter=Segment(splitter.end, splitter.start))
        other.validate()
        assert other.area == polygon.area - original_part.area
        LOGGER.debug("Second part of the polygon is updated.")

        if splitter_vertices[idx + 1] in border_vertices:
            border_vertices[border_vertices.index(
                splitter_vertices[idx + 1])] = splitter.start
        if splitter_vertices[idx - 1] in border_vertices:
            if splitter_vertices[idx - 1] not in polygon.border.vertices:
                border_vertices[border_vertices.index(
                    splitter_vertices[idx - 1])] = splitter.end
            else:
                if splitter.end != splitter_vertices[idx - 1]:
                    if orientation is Orientation.COUNTERCLOCKWISE:
                        border_vertices.insert(border_vertices.index(
                            splitter_vertices[idx - 1]) + 1, splitter.end)
                    else:
                        border_vertices.insert(border_vertices.index(
                            splitter_vertices[idx - 1]), splitter.end)

        splitter_vertices[idx - 1] = splitter.end
        splitter_vertices[idx + 1] = splitter.start
        if splitter.start == splitter_vertices[idx + 2]:
            splitter_vertices.pop(idx + 2)
        splitter_vertices.pop(idx)
        LOGGER.debug("Splitter vertices updated.")


def to_ccw_segments(*,
                    area: Fraction,
                    border_vertices: List[Point[Fraction]],
                    tail_start: Point[Fraction],
                    tail_vertex: Point[Fraction],
                    head_end: Point[Fraction],
                    head_vertex: Point[Fraction]) -> Tuple[Segment[Fraction],
                                                           Segment[Fraction],
                                                           Fraction]:
    """
    :param area:
    :param border_vertices:
    :param tail_start:
    :param tail_vertex:
    :param head_end:
    :param head_vertex:
    :return:
    """
    if (context.angle_orientation(tail_start, head_vertex, head_end)
            is Orientation.COLLINEAR):
        head_start = head_vertex
        old_tail_start = tail_start
        tail_start = to_ccw_end(area,
                                head_start,
                                tail_start,
                                tail_vertex)
        tail_end = to_ccw_end(area,
                              head_end,
                              old_tail_start,
                              tail_vertex)
        tail_segment = Segment(tail_start, tail_end)
        head_segment = Segment(head_start, head_end)
        area_diff = Polygon(Contour([head_end, old_tail_start, tail_start])
                            ).area
        return tail_segment, head_segment, area_diff
    head_start = to_cw_end(area,
                           tail_start,
                           head_end,
                           head_vertex)
    LOGGER.debug("Found head_start.")
    if tail_start not in border_vertices:
        LOGGER.debug("Tail_start is not in border vertices.")
        tail_end = to_ccw_end(area,
                              head_end,
                              tail_start,
                              tail_vertex)
    else:
        LOGGER.debug("Tail_start is on border vertices.")
        next_border_vertex = border_vertices[(border_vertices.index(
            tail_start) + 1) % len(border_vertices)]
        if context.angle_orientation(next_border_vertex,
                                     tail_start,
                                     head_end) is not Orientation.COLLINEAR:
            tail_end = to_ccw_end(area,
                                  head_end,
                                  tail_start,
                                  next_border_vertex)
        else:
            tail_end = next_border_vertex
            head_end = to_cw_end(area,
                                 tail_end,
                                 head_end,
                                 head_start)
    LOGGER.debug("Found tail_end.")
    tail_segment = Segment(tail_start, tail_end)
    head_segment = Segment(head_start, head_end)
    return tail_segment, head_segment, 0


def to_cw_segments(*,
                   area: Fraction,
                   tail_end: Point[Fraction],
                   tail_vertex: Point[Fraction],
                   head_start: Point[Fraction],
                   head_vertex: Point[Fraction]) -> Tuple[Segment[Fraction],
                                                          Segment[Fraction],
                                                          Fraction]:
    """
    :param area:
    :param tail_end:
    :param tail_vertex:
    :param head_start:
    :param head_vertex:
    :return:
    """
    if (context.angle_orientation(tail_end, head_start, head_vertex)
            is Orientation.COLLINEAR):
        LOGGER.debug("Head is collinear with triangle base.")
        head_end = head_vertex
        old_tail_end = tail_end
        tail_end = to_cw_end(area,
                             head_end,
                             tail_end,
                             tail_vertex)
        LOGGER.debug("New tail_end found.")
        tail_start = to_cw_end(area,
                               head_start,
                               old_tail_end,
                               tail_vertex)
        tail_segment = Segment(tail_start, tail_end)
        head_segment = Segment(head_start, head_end)
        area_diff = Polygon(Contour([head_start, old_tail_end, tail_end])).area
        return tail_segment, head_segment, area_diff
    tail_start = to_cw_end(area,
                           head_start,
                           tail_end,
                           tail_vertex)
    head_end = to_ccw_end(area,
                          tail_end,
                          head_start,
                          head_vertex)
    if (head_vertex in Segment(head_start, head_end)
            and head_vertex != head_end):
        head_end = head_vertex
        tail_end = to_cw_end(area - Polygon(Contour([head_start,
                                                     head_end,
                                                     tail_end])).area,
                             head_end,
                             tail_end,
                             tail_vertex)
    LOGGER.debug("Tail_start found.")
    tail_segment = Segment(tail_start, tail_end)
    head_segment = Segment(head_start, head_end)
    return tail_segment, head_segment, 0


def to_cw_end(area: Fraction,
              pivot: Point,
              origin: Point,
              destination: Point) -> Point:
    """
    :param area:
    :param pivot:
    :param origin:
    :param destination:
    :return:
    """
    if destination.x == origin.x:
        x_coordinate = origin.x
        y_coordinate = origin.y + 2 * area / (pivot.x - origin.x)
    else:
        slope = (destination.y - origin.y) / (destination.x - origin.x)
        intercept = origin.y - origin.x * slope
        x_coordinate = ((2 * area + pivot.x * origin.y - origin.x * pivot.y
                         - intercept * (pivot.x - origin.x))
                        / (origin.y - pivot.y + slope * (pivot.x - origin.x)))
        y_coordinate = slope * x_coordinate + intercept
    point = Point(x_coordinate, y_coordinate)
    assert Polygon(Contour([pivot, origin, point]), []).area == area
    return point


def to_ccw_end(area: Fraction,
               pivot: Point,
               origin: Point,
               destination: Point) -> Point:
    """
    :param area:
    :param pivot:
    :param origin:
    :param destination:
    :return:
    """
    if destination.x == origin.x:
        x_coordinate = origin.x
        y_coordinate = origin.y + 2 * area / (origin.x - pivot.x)
    else:
        slope = (destination.y - origin.y) / (destination.x - origin.x)
        intercept = origin.y - origin.x * slope
        x_coordinate = ((2 * area + origin.x * pivot.y - pivot.x * origin.y
                         - intercept * (origin.x - pivot.x))
                        / (pivot.y - origin.y + slope * (origin.x - pivot.x)))
        y_coordinate = slope * x_coordinate + intercept
    point = Point(x_coordinate, y_coordinate)
    assert Polygon(Contour([pivot, origin, point]), []).area == area
    return point


def update_part(part: Polygon,
                *,
                splitter_vertices: List[Point[Fraction]],
                idx: int,
                polygon_border_vertices: List[Point[Fraction]],
                border_vertices: List[Point[Fraction]],
                splitter: Segment[Fraction]):
    """
    :param part:
    :param splitter_vertices:
    :param idx:
    :param polygon_border_vertices:
    :param border_vertices:
    :param splitter:
    :return:
    """
    vertices = list(part.border.to_counterclockwise().vertices)
    i = vertices.index(splitter_vertices[idx])
    if vertices[i - 1] in polygon_border_vertices:
        j = border_vertices.index(vertices[i - 1])
        if (context.angle_orientation(
                border_vertices[j],
                border_vertices[(j + 1) % len(border_vertices)],
                vertices[i]) is Orientation.COLLINEAR
                and border_vertices[(j + 1) % len(border_vertices)]
                != vertices[i]):
            vertices.insert(i, border_vertices[(j + 1) % len(border_vertices)])
            vertices.insert(i + 1, splitter.end)
            i += 2
        elif (context.angle_orientation(
                splitter.end,
                border_vertices[j],
                border_vertices[(j + 1) % len(border_vertices)])
                is Orientation.COLLINEAR
                and splitter.end != border_vertices[j]):
            vertices.insert(i, splitter.end)
            i += 1
        else:
            vertices[i - 1] = splitter.end
    elif vertices[i - 2] == splitter.end:
        vertices.pop(i - 1)
        i -= 1
    else:
        vertices[i - 1] = splitter.end
    if vertices[(i + 1) % len(vertices)] in polygon_border_vertices:
        j = polygon_border_vertices.index(vertices[(i + 1) % len(vertices)])
        if (context.angle_orientation(vertices[i],
                                      polygon_border_vertices[j - 1],
                                      polygon_border_vertices[j])
                is Orientation.COLLINEAR):
            if splitter.start not in Segment(
                    polygon_border_vertices[j],
                    polygon_border_vertices[(j + 1)
                                            % len(polygon_border_vertices)]):
                vertices[i] = splitter.start
                vertices.insert(i + 1, polygon_border_vertices[j - 1])
            else:
                vertices[i] = splitter.start
                vertices.pop((i + 1) % len(vertices))
        elif (context.angle_orientation(splitter.start,
                                        vertices[(i + 1) % len(vertices)],
                                        vertices[(i + 2) % len(vertices)])
                is not Orientation.COLLINEAR):
            vertices[i] = splitter.start
        else:
            vertices[(i + 1) % len(vertices)] = splitter.start
            vertices.pop(i)
    elif vertices[(i + 2) % len(vertices)] == splitter.start:
        vertices.pop(i + 1)
        vertices.pop(i)
    else:
        vertices[(i + 1) % len(vertices)] = splitter.start
        vertices.pop(i)
    return Polygon(Contour(vertices), part.holes)
