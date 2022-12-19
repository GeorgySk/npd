import logging
import math
from fractions import Fraction

from gon.base import Polygon
from hypothesis import (Verbosity,
                        given,
                        settings)

from npd.draw import draw
from npd.npd import split_into_two
from tests.strategies.base import (polygons,
                                   steiner_points_counts,
                                   unit_area_fractions)

LOGGER = logging.getLogger(__name__)

draw.DRAW = False


@given(polygon=polygons,
       unit_area_fraction=unit_area_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_unions_equality(polygon: Polygon[Fraction],
                               unit_area_fraction: Fraction,
                               steiner_points_count: int) -> None:
    requirement = polygon.area * unit_area_fraction
    first_part, second_part = split_into_two(
        polygon,
        requirement=requirement,
        steiner_points_count=steiner_points_count)
    draw.DRAW = True
    draw(first_part, blue=second_part)
    draw.DRAW = False

    assert True
    # union = first_part | second_part
    # note(f"Resulting union: {union}\n")
    # assert polygon == first_part | second_part


@given(polygon=polygons,
       unit_area_fraction=unit_area_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_validity(polygon: Polygon[Fraction],
                        unit_area_fraction: Fraction,
                        steiner_points_count: int) -> None:
    if not hasattr(test_nodes_unions_equality, 'count'):
        test_nodes_unions_equality.count = 1
    else:
        test_nodes_unions_equality.count += 1
    requirement = polygon.area * unit_area_fraction
    draw(polygon)
    first_part, second_part = split_into_two(
        polygon,
        requirement=requirement,
        steiner_points_count=steiner_points_count)
    draw(first_part, blue=second_part)

    assert isinstance(first_part, Polygon)
    assert isinstance(second_part, Polygon)
