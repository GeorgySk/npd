from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.npd import to_partition
from tests.strategies.base import (sparse_deltas,
                                   varied_polygons)


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_triangles_sharing_an_edge(polygon: Polygon[Fraction],
                                   delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    for triangle in partition.triangular_map:
        result = partition.triangle_neighbors(triangle)

        expected = {other for other in partition.triangular_map
                    if other is not triangle
                    and set(triangle.edges) & set(other.edges)}
        assert len(result) == len(set(result)) <= len(triangle.edges)
        assert set(result) == expected


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_symmetry(polygon: Polygon[Fraction],
                  delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    assert all(triangle in partition.triangle_neighbors(neighbor)
               for triangle in partition.triangular_map
               for neighbor in partition.triangle_neighbors(triangle))
