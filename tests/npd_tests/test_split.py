import math
from fractions import Fraction
from typing import List

import matplotlib.pyplot as plt
from gon.base import Polygon
from hypothesis import (Verbosity,
                        given,
                        note,
                        settings)

from draw import draw
from npd.npd import (split,
                     unite)
from tests.strategies.base import (polygons,
                                   requirements_fractions,
                                   steiner_points_counts)

DRAW = True


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
@settings(verbosity=Verbosity.verbose)
def test_nodes_unions_equality(polygon: Polygon[Fraction],
                               requirements_fractions: List[Fraction],
                               steiner_points_count: Fraction) -> None:
    if not hasattr(test_nodes_unions_equality, 'count'):
        test_nodes_unions_equality.count = 1
    else:
        test_nodes_unions_equality.count += 1
    note(f"Polygon Nº{test_nodes_unions_equality.count}: {polygon}")
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    steiner_distance = Fraction(math.sqrt(polygon.area / steiner_points_count))
    parts = split(polygon,
                  requirements=requirements,
                  steiner_distance=steiner_distance)
    if DRAW:
        print("DRAWING")
        for part in parts:
            draw(part, fill=True, alpha=0.5)
        plt.title(str(test_nodes_unions_equality.count))
        print("FINISHED DRAWING")

    union = unite(parts)
    note(f"Resulting union: {union}")
    assert polygon == union


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
@settings(verbosity=Verbosity.verbose)
def test_nodes_validity(polygon: Polygon[Fraction],
                        requirements_fractions: List[Fraction],
                        steiner_points_count: Fraction) -> None:
    if not hasattr(test_nodes_unions_equality, 'count'):
        test_nodes_unions_equality.count = 1
    else:
        test_nodes_unions_equality.count += 1
    note(f"Polygon Nº{test_nodes_unions_equality.count}: {polygon}")
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    steiner_distance = Fraction(math.sqrt(polygon.area / steiner_points_count))
    if DRAW:
        draw(polygon)
    parts = split(polygon,
                  requirements=requirements,
                  steiner_distance=steiner_distance)
    if DRAW:
        for part in parts:
            draw(part, fill=True, alpha=0.5)
        plt.title(str(test_nodes_unions_equality.count))
    for part in parts:
        assert isinstance(part, Polygon)
