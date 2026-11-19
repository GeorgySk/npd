from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.npd import split
from npd.utils import unite
from tests.strategies.base import (polygons,
                                   requirements_fractions,
                                   steiner_points_counts)
from tests.utils import (HOLED_SQUARE,
                         NON_CONVEX_POLYGONS)


def assert_valid_split(polygon: Polygon[Fraction],
                       requirements: list[Fraction],
                       parts: list[Polygon[Fraction]]) -> None:
    assert all(isinstance(part, Polygon) for part in parts)
    # `split` cuts off the requirements in increasing order.
    assert [part.area for part in parts] == sorted(requirements)
    # With the areas adding up, an equal union means the parts don't
    # overlap.
    assert unite(parts) == polygon


# Known failure: the digits limit for integer string conversion (see
# CLAUDE.md).
@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
def test_split(polygon: Polygon[Fraction],
               requirements_fractions: list[Fraction],
               steiner_points_count: int) -> None:
    requirements = [fraction * polygon.area
                    for fraction in requirements_fractions]

    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)

    assert_valid_split(polygon, requirements, parts)


@pytest.mark.parametrize("polygon", NON_CONVEX_POLYGONS)
@pytest.mark.parametrize("requirements_fractions",
                         [[Fraction(1, 3)] * 3,
                          [Fraction(1, 2), Fraction(1, 6), Fraction(1, 3)],
                          [Fraction(n, 10) for n in range(1, 5)]])
@pytest.mark.parametrize("steiner_points_count", [4, 16])
def test_non_convex(polygon: Polygon[Fraction],
                    requirements_fractions: list[Fraction],
                    steiner_points_count: int) -> None:
    requirements = [fraction * polygon.area
                    for fraction in requirements_fractions]

    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)

    assert_valid_split(polygon, requirements, parts)


def test_remainder_with_hole_touching_border() -> None:
    # The second split's remainder has its hole touching the outer border
    # at (3, 3/2); its triangulation is empty.
    requirements = [HOLED_SQUARE.area * Fraction(n, 10) for n in range(1, 5)]

    parts = split(HOLED_SQUARE,
                  requirements=requirements,
                  steiner_points_count=4)

    assert_valid_split(HOLED_SQUARE, requirements, parts)


@given(polygon=polygons,
       steiner_points_count=steiner_points_counts)
def test_single_requirement(polygon: Polygon[Fraction],
                            steiner_points_count: int) -> None:
    assert split(polygon,
                 requirements=[polygon.area],
                 steiner_points_count=steiner_points_count) == [polygon]


@given(polygon=polygons,
       requirements_fractions=requirements_fractions,
       steiner_points_count=steiner_points_counts)
def test_requirements_not_summing_to_area(
        polygon: Polygon[Fraction],
        requirements_fractions: list[Fraction],
        steiner_points_count: int) -> None:
    requirements = [fraction * polygon.area
                    for fraction in requirements_fractions[1:]]

    with pytest.raises(ValueError, match="equal to polygon area"):
        split(polygon,
              requirements=requirements,
              steiner_points_count=steiner_points_count)
