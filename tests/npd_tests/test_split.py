from fractions import Fraction
from typing import List

import matplotlib.pyplot as plt
from gon.base import Polygon
from hypothesis import (Verbosity,
                        given,
                        settings)

from npd.draw import draw
from npd.npd import (split,
                     unite)
from tests.strategies.base import (polygons,
                                   requirements_fractions,
                                   steiner_points_counts)

draw.DRAW = False


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_unions_equality(polygon: Polygon[Fraction],
                               requirements_fractions: List[Fraction],
                               steiner_points_count: int) -> None:
    if not hasattr(test_nodes_unions_equality, 'count'):
        test_nodes_unions_equality.count = 1
    else:
        test_nodes_unions_equality.count += 1
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)
    if draw.DRAW:
        for part in parts:
            draw(part, fill=True, alpha=0.5)
        plt.title(str(test_nodes_unions_equality.count))

    union = unite(parts)
    assert polygon == union


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
def test_nodes_validity(polygon: Polygon[Fraction],
                        requirements_fractions: List[Fraction],
                        steiner_points_count: int) -> None:
    if not hasattr(test_nodes_unions_equality, 'count'):
        test_nodes_unions_equality.count = 1
    else:
        test_nodes_unions_equality.count += 1
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    draw(polygon)
    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)
    if draw.DRAW:
        for part in parts:
            draw(part, fill=True, alpha=0.5)
        plt.title(str(test_nodes_unions_equality.count))
    for part in parts:
        assert isinstance(part, Polygon)
