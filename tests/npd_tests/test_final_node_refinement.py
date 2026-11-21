import math
from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import final_node_refinement
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition,
                            to_triangular_map)
from tests.strategies.base import (partition_operations,
                                   refinements,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         border_length,
                         is_connected,
                         rebuilt_neighbors)

Setup = tuple[Chunk, Chunk, CachedPolygon]


def to_partition_after(polygon: Polygon[Fraction],
                       delta: Fraction,
                       operations: list[tuple[bool, int, int]]) -> Partition:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    return partition


def setups(partition: Partition) -> list[Setup]:
    """
    Lists the cases ``readjust`` refines.

    A node of the chunk touching the neighbor, with a neighbor in the
    chunk, whose removal leaves the chunk connected.
    """
    return [(chunk, neighbor, node)
            for chunk in partition.chunk_map
            for node in chunk.triangles
            if any(other in chunk.triangles
                   for other in partition.triangle_neighbors(node))
            and not partition.is_articulation_triangle(node, chunk.triangles)
            for neighbor in partition.chunk_neighbors(chunk)
            if any(other in neighbor.triangles
                   for other in partition.triangle_neighbors(node))]


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       refinement=refinements)
def test_meets_requirement(polygon: Polygon[Fraction],
                           delta: Fraction,
                           operations: list[tuple[bool, int, int]],
                           refinement: tuple[int, Fraction]) -> None:
    partition = to_partition_after(polygon, delta, operations)
    choice, fraction = refinement
    cases = setups(partition)
    if not cases:
        return
    chunk, neighbor, node = cases[choice % len(cases)]
    requirement = chunk.area - node.area + fraction * node.area
    total = chunk.area + neighbor.area
    connected = [chunk_ for chunk_ in (chunk, neighbor)
                 if is_connected(partition, chunk_.triangles)]

    final_node_refinement(node,
                          requirement=requirement,
                          chunk=chunk,
                          neighbor=neighbor,
                          partition=partition)

    assert chunk.area == requirement
    assert neighbor.area == total - requirement
    assert all(is_connected(partition, chunk_.triangles)
               for chunk_ in connected)
    assert partition.chunk_map == rebuilt_neighbors(partition)
    assert partition.triangular_map == to_triangular_map(
        list(partition.triangular_map))
    for chunk_ in partition.chunk_map:
        assert chunk_.area == sum(triangle.area
                                  for triangle in chunk_.triangles)
        assert math.isclose(chunk_.perimeter,
                            border_length(chunk_.triangles, partition))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       refinement=refinements)
def test_rejects_unreachable_requirement(polygon: Polygon[Fraction],
                                         delta: Fraction,
                                         operations: list[tuple[bool,
                                                                int,
                                                                int]],
                                         refinement: tuple[int, Fraction]
                                         ) -> None:
    partition = to_partition_after(polygon, delta, operations)
    choice, _ = refinement
    cases = setups(partition)
    if not cases:
        return
    chunk, neighbor, node = cases[choice % len(cases)]

    for requirement in (chunk.area, chunk.area - node.area):
        with pytest.raises(ValueError, match="requirement"):
            final_node_refinement(node,
                                  requirement=requirement,
                                  chunk=chunk,
                                  neighbor=neighbor,
                                  partition=partition)
