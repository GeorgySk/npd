import logging
import math
from fractions import Fraction
from itertools import (chain,
                       product)
from typing import (List,
                    Tuple)

from gon.base import (Point,
                      Polygon)
from ground.base import get_context

from npd import (convex,
                 nonconvex)
from npd.draw import draw
from npd.structures import Partition
from npd.utils import (fast_length,
                       unite)

LOGGER = logging.getLogger(__name__)
context = get_context()


def split(polygon: Polygon[Fraction],
          *,
          requirements: List[Fraction],
          steiner_points_count: int) -> List[Polygon[Fraction]]:
    """
    Splits a polygon into N parts where N is length of the list of
    area requirements.
    :param polygon: input polygon, with or without holes
    :param requirements: list of area requirements (its sum must be
    equal to the area of the input polygon)
    :param steiner_points_count: approximate number of steiner points
    that will be used for splitting. Larger values will result in less
    time to compute but worse results in terms of compactness of the
    resulting parts
    :return: a list of resulting polygon parts
    """
    if not sum(requirements) == polygon.area:
        raise ValueError("Sum of requirements must be equal to polygon area")
    result = []
    remainder = polygon
    for requirement in sorted(requirements)[:-1]:
        part, remainder = split_into_two(
            remainder,
            requirement=requirement,
            steiner_points_count=steiner_points_count)
        result.append(part)
        draw(*result, remainder)
    result.append(remainder)
    return result


def split_into_two(polygon: Polygon[Fraction],
                   *,
                   requirement: Fraction,
                   steiner_points_count: int
                   ) -> Tuple[Polygon[Fraction], Polygon[Fraction]]:
    """
    Splits a polygon into two parts according to the given area
    requirement
    :param polygon: input polygon
    :param requirement: input area requirement (should be less or
    equal than the half of the polygon's area)
    :param steiner_points_count: an approximate number of steiner points
    that will be used for splitting. Larger values will result in less
    time to compute but worse results in terms of compactness of the
    resulting parts
    :return: a polygon part with the specified area and the remaining
    part
    """
    if requirement > Fraction(1, 2) * polygon.area:
        raise ValueError("Requirement can't be larger than half of area.")
    if polygon.is_convex:
        part, other = convex.split(contour=polygon.border,
                                   area_requirement=requirement,
                                   key=lambda x, y: x.length)
        return Polygon(part, []), Polygon(other, [])
    delta = Fraction(math.sqrt(polygon.area / steiner_points_count))
    partition = to_partition(polygon, delta=delta)
    draw(*partition.triangular_map)
    draw(*partition.chunk_map)
    chunk, remainder = nonconvex.split(partition, requirement=requirement)
    part = unite(chunk.triangles)
    other = unite(remainder.triangles)
    return part, other


def to_partition(polygon: Polygon,
                 *,
                 delta: Fraction) -> Partition:
    """
    Converts the given polygon to a Partition object
    :param polygon: input polygon
    :param delta: an approximate distance between Steiner points to be
    inserted in the polygon when performing triangulation
    :return: Partition object
    """
    extra_points = steiner_points(polygon, delta=delta)
    draw(polygon, *extra_points)
    return Partition(polygon, extra_points=extra_points)


def steiner_points(polygon: Polygon[Fraction],
                   delta: Fraction) -> List[Point[Fraction]]:
    """
    :param polygon:
    :param delta:
    :return:
    """
    box = context.polygon_box(polygon)
    xs_grid_count = (box.max_x - box.min_x) // delta
    ys_grid_count = (box.max_y - box.min_y) // delta
    offsets_product = product(range(1, xs_grid_count), range(1, ys_grid_count))
    grid = [Point(box.min_x + x * delta, box.min_y + y * delta)
            for x, y in offsets_product]
    constrained_grid = [point for point in grid if point in polygon]
    points_on_border = []
    for segment in chain(polygon.border.segments,
                         *[hole.segments for hole in polygon.holes]):
        extra_points_count = Fraction(fast_length(segment)) // delta
        if not extra_points_count:
            continue
        delta_x = (segment.end.x - segment.start.x) / extra_points_count
        delta_y = (segment.end.y - segment.start.y) / extra_points_count
        for i in range(1, extra_points_count + 1):
            points_on_border.append(Point(segment.start.x + i * delta_x,
                                          segment.start.y + i * delta_y))
    return constrained_grid + points_on_border
