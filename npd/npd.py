"""Splitting of polygons into parts with the given areas."""
import math
from fractions import Fraction
from itertools import (chain,
                       product)

from gon.base import (Point,
                      Polygon)
from ground.base import get_context
from pdan import minimizing_split

from npd import nonconvex
from npd.structures import Partition
from npd.utils import (fast_length,
                       unite)

context = get_context()


def split(polygon: Polygon[Fraction],
          *,
          requirements: list[Fraction],
          steiner_points_count: int) -> list[Polygon[Fraction]]:
    """
    Splits a polygon into N parts where N is the number of requirements.

    Cuts off one part at a time with ``split_into_two``, smallest
    requirement first: the smallest of the requirements left is at most
    half of the remainder's area, as ``split_into_two`` needs. The last
    remainder is the part for the largest requirement.
    :param polygon: input polygon, with or without holes
    :param requirements: positive area requirements adding up to the
    polygon's area
    :param steiner_points_count: approximate number of Steiner points to
    triangulate each non-convex polygon with (see ``split_into_two``);
    larger values take more time and give more compact parts
    :return: the parts, the i-th with the i-th requirement's area
    :raises ValueError: if the requirements don't add up to the area or
    aren't positive
    """
    if sum(requirements) != polygon.area:
        raise ValueError("Sum of requirements must be equal to polygon area")
    order = sorted(range(len(requirements)), key=requirements.__getitem__)
    parts: list[Polygon[Fraction]] = [polygon] * len(requirements)
    remainder = polygon
    for index in order[:-1]:
        parts[index], remainder = split_into_two(
            remainder,
            requirement=requirements[index],
            steiner_points_count=steiner_points_count)
    parts[order[-1]] = remainder
    return parts


def split_into_two(polygon: Polygon[Fraction],
                   *,
                   requirement: Fraction,
                   steiner_points_count: int
                   ) -> tuple[Polygon[Fraction], Polygon[Fraction]]:
    """
    Splits a polygon into two parts according to the area requirement.

    A convex polygon is cut by a single segment with ``pdan``, keeping
    the part with the least perimeter; any other polygon is triangulated
    with Steiner points spaced ``sqrt(area / steiner_points_count)``
    apart and split by ``nonconvex.split``. Both parts come out as
    single polygons: ``nonconvex.split``'s parts are connected through
    shared edges, so their triangles unite into one polygon each.
    :param polygon: input polygon
    :param requirement: area of the part, positive and at most half of
    the polygon's area
    :param steiner_points_count: approximate number of Steiner points to
    triangulate a non-convex polygon with; larger values take more time
    and give more compact parts
    :return: a polygon part with the specified area and the remaining
    part
    :raises ValueError: if the requirement isn't positive and at most
    half of the area, or the Steiner points count isn't positive
    """
    if not 0 < requirement <= Fraction(1, 2) * polygon.area:
        raise ValueError("Requirement must be positive and at most half "
                         "of the area.")
    if steiner_points_count <= 0:
        raise ValueError("Steiner points count must be positive.")
    if polygon.is_convex:
        # gon's exact lengths are symbolic: slower, and comparing them
        # converts ever larger integers to strings, which hits Python's
        # digits limit over repeated splits.
        part, other = minimizing_split(
            contour=polygon.border,
            area_requirement=requirement,
            key=lambda contour, _: fast_length(contour))
        return Polygon(part), Polygon(other)
    delta = Fraction(math.sqrt(polygon.area / steiner_points_count))
    partition = to_partition(polygon, delta=delta)
    chunk, remainder = nonconvex.split(partition, requirement=requirement)
    part = unite(chunk.triangles)
    other = unite(remainder.triangles)
    return part, other


def to_partition(polygon: Polygon[Fraction],
                 *,
                 delta: Fraction) -> Partition:
    """
    Converts the given polygon to a Partition object.

    :param polygon: input polygon
    :param delta: an approximate distance between Steiner points to be
    inserted in the polygon when performing triangulation
    :return: Partition object
    """
    extra_points = steiner_points(polygon, delta=delta)
    return Partition(polygon, extra_points=extra_points)


def steiner_points(polygon: Polygon[Fraction],
                   delta: Fraction) -> list[Point[Fraction]]:
    """
    Generates Steiner points inside and on the border of the polygon.

    :param polygon: input polygon
    :param delta: an approximate distance between the points
    :return: distinct points that aren't vertices of the polygon
    """
    box = context.polygon_box(polygon)
    xs_grid_count = (box.max_x - box.min_x) // delta
    ys_grid_count = (box.max_y - box.min_y) // delta
    offsets_product = product(range(1, xs_grid_count), range(1, ys_grid_count))
    grid = [Point(box.min_x + x * delta, box.min_y + y * delta)
            for x, y in offsets_product]
    constrained_grid = [point for point in grid if point in polygon]
    points_on_border: list[Point[Fraction]] = []
    for segment in chain(polygon.border.segments,
                         *[hole.segments for hole in polygon.holes]):
        extra_points_count = Fraction(fast_length(segment)) // delta
        if not extra_points_count:
            continue
        delta_x = (segment.end.x - segment.start.x) / extra_points_count
        delta_y = (segment.end.y - segment.start.y) / extra_points_count
        points_on_border.extend(Point(segment.start.x + i * delta_x,
                                      segment.start.y + i * delta_y)
                                for i in range(1, extra_points_count))
    vertices = set(chain(polygon.border.vertices,
                         *[hole.vertices for hole in polygon.holes]))
    # Grid points can land on the border, repeating a border point or
    # hitting a polygon vertex; `dict.fromkeys` drops repeats in order.
    return [point
            for point in dict.fromkeys(chain(constrained_grid,
                                             points_on_border))
            if point not in vertices]
