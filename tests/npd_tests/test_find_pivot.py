from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon,
                      Segment)
from hypothesis import given

from npd.nonconvex import find_pivot
from npd.npd import to_partition
from npd.structures import Partition
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import apply


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_pivot_edges(polygon: Polygon[Fraction],
                     delta: Fraction,
                     operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        for triangle, neighbors in partition.triangular_map.items():
            result = find_pivot(triangle, partition=partition, chunk=chunk)

            vertices = triangle.border.vertices
            before = Segment(vertices[result - 1], vertices[result])
            after = Segment(vertices[result], vertices[result - 2])
            across = [neighbor for neighbor in neighbors.values()
                      if neighbor is not None]
            inside = [neighbor for neighbor in across
                      if neighbor in chunk.triangles]
            if inside:
                assert neighbors[before] is inside[0]
            elif across:
                assert neighbors[after] is across[0]
            else:
                assert result == 0


def test_lone_triangle() -> None:
    border = Contour([Point(Fraction(0), Fraction(0)),
                      Point(Fraction(1), Fraction(0)),
                      Point(Fraction(0), Fraction(1))])
    partition = Partition(Polygon(border, []))
    (chunk,) = partition.chunk_map
    (triangle,) = chunk.triangles

    assert find_pivot(triangle, partition=partition, chunk=chunk) == 0
