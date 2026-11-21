import math
from fractions import Fraction

import pytest
from gon.base import (Point,
                      Polygon)
from hypothesis import given

from npd.nonconvex import divide
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Partition,
                            to_triangular_map)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   triangle_cuts,
                                   varied_polygons)
from tests.utils import (apply,
                         border_length,
                         rebuilt_neighbors)


def to_new_vertex(node: CachedPolygon,
                  pivot_index: int,
                  fraction: Fraction) -> Point[Fraction]:
    """Point at the given fraction of the edge opposite the pivot."""
    vertices = node.border.vertices
    start, end = vertices[pivot_index - 1], vertices[pivot_index - 2]
    return Point(start.x + fraction * (end.x - start.x),
                 start.y + fraction * (end.y - start.y))


def to_partition_after(polygon: Polygon[Fraction],
                       delta: Fraction,
                       operations: list[tuple[bool, int, int]]) -> Partition:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    return partition


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       cut=triangle_cuts)
def test_keeps_partition_consistent(polygon: Polygon[Fraction],
                                   delta: Fraction,
                                   operations: list[tuple[bool, int, int]],
                                   cut: tuple[int, int, Fraction]) -> None:
    partition = to_partition_after(polygon, delta, operations)
    choice, pivot_index, fraction = cut
    triangles = list(partition.triangular_map)
    node = triangles[choice % len(triangles)]
    chunk = next(chunk for chunk in partition.chunk_map
                 if node in chunk.triangles)
    areas = {chunk_: chunk_.area for chunk_ in partition.chunk_map}

    left, right = divide(partition,
                         node=node,
                         pivot_index=pivot_index,
                         new_vertex=to_new_vertex(node, pivot_index,
                                                  fraction),
                         chunk=chunk)

    assert left.area + right.area == node.area
    assert node not in partition.triangular_map
    assert left in chunk.triangles
    assert right in chunk.triangles
    assert partition.triangular_map == to_triangular_map(
        list(partition.triangular_map))
    assert partition.chunk_map == rebuilt_neighbors(partition)
    owned = [triangle
             for chunk_ in partition.chunk_map
             for triangle in chunk_.triangles]
    assert len(owned) == len(partition.triangular_map)
    assert set(owned) == set(partition.triangular_map)
    for chunk_ in partition.chunk_map:
        assert chunk_.area == areas[chunk_]
        assert chunk_.area == sum(triangle.area
                                  for triangle in chunk_.triangles)
        assert math.isclose(chunk_.perimeter,
                            border_length(chunk_.triangles, partition))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_rejects_node_of_other_chunk(polygon: Polygon[Fraction],
                                     delta: Fraction,
                                     operations: list[tuple[bool, int, int]]
                                     ) -> None:
    partition = to_partition_after(polygon, delta, operations)
    chunk, *rest = partition.chunk_map
    if not rest:
        return
    node = next(iter(rest[0].triangles))

    with pytest.raises(ValueError, match="belong to the chunk"):
        divide(partition,
               node=node,
               pivot_index=0,
               new_vertex=to_new_vertex(node, 0, Fraction(1, 2)),
               chunk=chunk)
