from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.npd import (steiner_points,
                     to_partition)
from tests.strategies.base import (sparse_deltas,
                                   varied_polygons)


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_steiner_points_used(polygon: Polygon[Fraction],
                             delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    vertices = {vertex
                for triangle in partition.triangular_map
                for vertex in triangle.border.vertices}
    assert set(steiner_points(polygon, delta)) <= vertices
