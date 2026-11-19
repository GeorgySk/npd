from fractions import Fraction

import pytest
from gon.base import (Contour,
                      Point,
                      Polygon)
from hypothesis import given

from npd.npd import split_into_two
from tests.strategies.base import (above_half_requirement_factors,
                                   half_requirement_factors,
                                   polygons,
                                   steiner_points_counts,
                                   varied_polygons)


def to_polygon(vertices: list[tuple[int, int]],
               *holes: list[tuple[int, int]]) -> Polygon[Fraction]:
    def to_contour(points: list[tuple[int, int]]) -> Contour[Fraction]:
        return Contour([Point(Fraction(x), Fraction(y)) for x, y in points])

    return Polygon(to_contour(vertices), [to_contour(hole) for hole in holes])


L_SHAPE = to_polygon([(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)])
U_SHAPE = to_polygon([(0, 0), (3, 0), (3, 2), (2, 2),
                      (2, 1), (1, 1), (1, 2), (0, 2)])
HOLED_SQUARE = to_polygon([(0, 0), (3, 0), (3, 3), (0, 3)],
                          [(1, 1), (1, 2), (2, 2), (2, 1)])


def assert_valid_split(polygon: Polygon[Fraction],
                       requirement: Fraction,
                       part: Polygon[Fraction],
                       other: Polygon[Fraction]) -> None:
    assert isinstance(part, Polygon)
    assert isinstance(other, Polygon)
    assert part.area == requirement
    assert other.area == polygon.area - requirement
    # With the areas adding up, an equal union means the parts don't
    # overlap.
    assert part | other == polygon


@given(polygon=varied_polygons,
       factor=half_requirement_factors,
       steiner_points_count=steiner_points_counts)
def test_split(polygon: Polygon[Fraction],
               factor: Fraction,
               steiner_points_count: int) -> None:
    requirement = polygon.area * factor

    part, other = split_into_two(polygon,
                                 requirement=requirement,
                                 steiner_points_count=steiner_points_count)

    assert_valid_split(polygon, requirement, part, other)


@pytest.mark.parametrize("polygon", [L_SHAPE, U_SHAPE, HOLED_SQUARE])
@pytest.mark.parametrize("factor", [Fraction(1, 10),
                                    Fraction(1, 3),
                                    Fraction(1, 2)])
@pytest.mark.parametrize("steiner_points_count", [4, 16])
def test_non_convex(polygon: Polygon[Fraction],
                    factor: Fraction,
                    steiner_points_count: int) -> None:
    requirement = polygon.area * factor

    part, other = split_into_two(polygon,
                                 requirement=requirement,
                                 steiner_points_count=steiner_points_count)

    assert_valid_split(polygon, requirement, part, other)


@given(polygon=polygons,
       factor=above_half_requirement_factors,
       steiner_points_count=steiner_points_counts)
def test_requirement_above_half(polygon: Polygon[Fraction],
                                factor: Fraction,
                                steiner_points_count: int) -> None:
    with pytest.raises(ValueError, match="half of the area"):
        split_into_two(polygon,
                       requirement=polygon.area * factor,
                       steiner_points_count=steiner_points_count)
