from fractions import Fraction

from gon.base import (EMPTY,
                      Multipolygon,
                      Polygon)
from hypothesis import given

from npd.npd import to_partition
from npd.utils import unite
from tests.strategies.base import (polygons,
                                   polygons_lists,
                                   polygons_lists_with_permutations,
                                   sparse_deltas)


def test_empty() -> None:
    assert unite([]) is EMPTY


@given(polygons)
def test_single(polygon: Polygon[Fraction]) -> None:
    assert unite([polygon]) == polygon


@given(polygons_lists)
def test_result_type(polygons_list: list[Polygon[Fraction]]) -> None:
    assert isinstance(unite(polygons_list), (Polygon, Multipolygon))


@given(polygons_lists)
def test_containment(polygons_list: list[Polygon[Fraction]]) -> None:
    result = unite(polygons_list)

    assert all(polygon <= result for polygon in polygons_list)


@given(polygons_lists)
def test_area_bounds(polygons_list: list[Polygon[Fraction]]) -> None:
    areas = [polygon.area for polygon in polygons_list]

    assert max(areas) <= unite(polygons_list).area <= sum(areas)


@given(polygons_lists_with_permutations)
def test_order_independence(
        polygons_list_with_permutation: tuple[list[Polygon[Fraction]],
                                              list[Polygon[Fraction]]]
) -> None:
    polygons_list, permutation = polygons_list_with_permutation

    assert unite(permutation) == unite(polygons_list)


@given(polygons)
def test_triangulation_union(polygon: Polygon[Fraction]) -> None:
    triangles = [Polygon(contour)
                 for contour in polygon.triangulate().triangles()]

    assert unite(triangles) == polygon


@given(polygons, sparse_deltas)
def test_partition_union(polygon: Polygon[Fraction],
                         delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    assert unite(partition.triangular_map) == polygon
