from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.structures import (CachedPolygon,
                            to_triangular_map)
from tests.strategies.base import (polygons_with_triangles,
                                   triangles_lists_with_permutations)


@given(polygons_with_triangles)
def test_keys(polygon_with_triangles: tuple[Polygon[Fraction],
                                            list[CachedPolygon]]) -> None:
    _, triangles = polygon_with_triangles

    result = to_triangular_map(triangles)

    assert list(result) == triangles
    assert all(set(result[triangle]) == set(triangle.edges)
               for triangle in triangles)


@given(polygons_with_triangles)
def test_neighbors_share_edge(
        polygon_with_triangles: tuple[Polygon[Fraction],
                                      list[CachedPolygon]]) -> None:
    _, triangles = polygon_with_triangles

    result = to_triangular_map(triangles)

    assert all(neighbor is None
               or (neighbor != triangle and edge in neighbor.edges)
               for triangle, neighbors in result.items()
               for edge, neighbor in neighbors.items())


@given(polygons_with_triangles)
def test_symmetry(polygon_with_triangles: tuple[Polygon[Fraction],
                                                list[CachedPolygon]]) -> None:
    _, triangles = polygon_with_triangles

    result = to_triangular_map(triangles)

    assert all(result[neighbor][edge] == triangle
               for triangle, neighbors in result.items()
               for edge, neighbor in neighbors.items()
               if neighbor is not None)


@given(polygons_with_triangles)
def test_unpaired_edges_on_boundary(
        polygon_with_triangles: tuple[Polygon[Fraction],
                                      list[CachedPolygon]]) -> None:
    polygon, triangles = polygon_with_triangles
    contours = [polygon.border, *polygon.holes]

    result = to_triangular_map(triangles)

    assert all((neighbor is None)
               is any(edge <= contour for contour in contours)
               for neighbors in result.values()
               for edge, neighbor in neighbors.items())


@given(triangles_lists_with_permutations)
def test_order_independence(
        triangles_with_permutation: tuple[list[CachedPolygon],
                                          list[CachedPolygon]]) -> None:
    triangles, permutation = triangles_with_permutation

    assert to_triangular_map(permutation) == to_triangular_map(triangles)
