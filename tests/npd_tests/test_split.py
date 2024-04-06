from fractions import Fraction
from typing import List

from gon.base import Polygon
from hypothesis import given

from npd.npd import (split,
                     unite)
from tests.strategies.base import (polygons,
                                   requirements_fractions,
                                   steiner_points_counts)


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_unions_equality(polygon: Polygon[Fraction],
                               requirements_fractions: List[Fraction],
                               steiner_points_count: int) -> None:
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)
    union = unite(parts)
    assert polygon == union


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_validity(polygon: Polygon[Fraction],
                        requirements_fractions: List[Fraction],
                        steiner_points_count: int) -> None:
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)
    for part in parts:
        assert isinstance(part, Polygon)
