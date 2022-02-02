from fractions import Fraction

from gon.base import Polygon
from hypothesis import (Verbosity,
                        given,
                        settings)

from npd.npd import to_partition
from tests.strategies.base import (polygons,
                                   steiner_points_counts)

DRAW = False


@given(polygon=polygons,
       steiner_points_count=steiner_points_counts)
@settings(verbosity=Verbosity.verbose)
def test_non_articulation_node_presence(polygon: Polygon[Fraction],
                                        steiner_points_count: int) -> None:
    partition = to_partition(polygon,
                             steiner_points_count=steiner_points_count)
    # There is no need to split partition into chunks and using only a
    # particular chunk for the test.
    assert any(
        not partition.is_articulation_triangle(
            triangle, chunk_triangles=partition.triangular_map)
        for triangle in partition.triangular_map)
