from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.npd import split_into_two
from tests.strategies.base import (polygons,
                                   steiner_points_counts,
                                   unit_area_fractions)


@given(polygon=polygons,
       unit_area_fraction=unit_area_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_unions_equality(polygon: Polygon[Fraction],
                               unit_area_fraction: Fraction,
                               steiner_points_count: int) -> None:
    requirement = polygon.area * unit_area_fraction
    split_into_two(polygon,
                   requirement=requirement,
                   steiner_points_count=steiner_points_count)


@given(polygon=polygons,
       unit_area_fraction=unit_area_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_validity(polygon: Polygon[Fraction],
                        unit_area_fraction: Fraction,
                        steiner_points_count: int) -> None:
    requirement = polygon.area * unit_area_fraction
    first_part, second_part = split_into_two(
        polygon,
        requirement=requirement,
        steiner_points_count=steiner_points_count)

    assert isinstance(first_part, Polygon)
    assert isinstance(second_part, Polygon)
