from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.npd import (steiner_points,
                     to_partition)
from npd.utils import unite
from tests.strategies.base import (sparse_deltas,
                                   varied_polygons)
from tests.utils import TOUCHING_HOLES_POLYGONS


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_steiner_points_used(polygon: Polygon[Fraction],
                             delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    vertices = {vertex
                for triangle in partition.triangular_map
                for vertex in triangle.border.vertices}
    assert set(steiner_points(polygon, delta)) <= vertices


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_union(polygon: Polygon[Fraction],
               delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    assert unite(partition.triangular_map) == polygon


@pytest.mark.parametrize("polygon",
                         [polygon for polygon, _ in TOUCHING_HOLES_POLYGONS])
@pytest.mark.parametrize("delta", [Fraction(1), Fraction(1, 3)])
def test_contours_touching_inside_edge(polygon: Polygon[Fraction],
                                       delta: Fraction) -> None:
    # gon's constrained Delaunay triangulation returns no triangles or
    # raises unless the touching point is passed as an extra point.
    partition = to_partition(polygon, delta=delta)

    assert unite(partition.triangular_map) == polygon
