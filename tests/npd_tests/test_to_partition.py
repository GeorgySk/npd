from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.npd import (steiner_points,
                     to_partition)
from npd.utils import unite
from tests.strategies.base import (sparse_deltas,
                                   varied_polygons)
from tests.utils import to_polygon


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_steiner_points_used(polygon: Polygon[Fraction],
                             delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    vertices = {vertex
                for triangle in partition.triangular_map
                for vertex in triangle.border.vertices}
    assert set(steiner_points(polygon, delta)) <= vertices


def test_hole_touching_border() -> None:
    # A hole may touch the border at a point; gon's constrained Delaunay
    # triangulation returns no triangles for such a polygon.
    polygon = to_polygon([(0, 0), (3, 0), (3, 3), (0, 3)],
                         [(1, 1), (1, 2), (3, Fraction(3, 2))])

    partition = to_partition(polygon, delta=Fraction(1))

    assert unite(partition.triangular_map) == polygon
