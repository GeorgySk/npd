from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon)
from hypothesis import given

from npd.utils import inverse_shoelace
from tests.strategies.base import (inner_fractions,
                                   triangle_vertices)

Vertices = tuple[Point[Fraction], Point[Fraction], Point[Fraction]]


def triangle_area(vertices: Vertices) -> Fraction:
    area: Fraction = Polygon(Contour(list(vertices)), []).area
    return area


@given(vertices=triangle_vertices,
       fraction=inner_fractions)
def test_area(vertices: Vertices, fraction: Fraction) -> None:
    endpoint, base_start, base_end = vertices
    area = fraction * triangle_area(vertices)

    result = inverse_shoelace(area,
                              endpoint=endpoint,
                              base_start=base_start,
                              base_end=base_end)

    assert triangle_area((endpoint, base_start, result)) == area


@given(vertices=triangle_vertices,
       fraction=inner_fractions)
def test_point_inside_base(vertices: Vertices, fraction: Fraction) -> None:
    endpoint, base_start, base_end = vertices
    area = fraction * triangle_area(vertices)

    result = inverse_shoelace(area,
                              endpoint=endpoint,
                              base_start=base_start,
                              base_end=base_end)

    # the area is that fraction of the whole triangle's, so the point
    # divides the base in that ratio
    assert result == Point(
        base_start.x + fraction * (base_end.x - base_start.x),
        base_start.y + fraction * (base_end.y - base_start.y))
