from fractions import Fraction

from gon.base import Polygon
from hypothesis import (Verbosity,
                        given,
                        settings)

from npd.npd import to_partition
from tests.strategies.base import (deltas,
                                   polygons)


@given(polygon=polygons,
       delta=deltas)
def test_non_articulation_node_presence(polygon: Polygon[Fraction],
                                        delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)
    # There is no need to split partition into chunks and using only a
    # particular chunk for the test.
    assert any(
        not partition.is_articulation_triangle(
            triangle, chunk_triangles=partition.triangular_map)
        for triangle in partition.triangular_map)
