from fractions import Fraction

import pytest
from gon.base import (Polygon,
                      Segment)
from hypothesis import given

from npd.npd import to_partition
from npd.structures import CachedPolygon
from tests.strategies.base import (sparse_deltas,
                                   varied_polygons)
from tests.utils import (to_inner_triangle_partition,
                         to_point)


def to_edges(triangle: CachedPolygon) -> set[Segment[Fraction]]:
    vertices = triangle.border.vertices
    return {Segment(vertices[index - 1], vertices[index])
            for index in range(len(vertices))}


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_shared_edge(polygon: Polygon[Fraction],
                     delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    for triangle in partition.triangular_map:
        for other in partition.triangle_neighbors(triangle):
            edge = partition.edge_between_triangles(triangle, other)

            assert to_edges(triangle) & to_edges(other) == {edge}
            assert partition.edge_between_triangles(other, triangle) == edge


@given(polygon=varied_polygons,
       delta=sparse_deltas)
def test_not_neighbors(polygon: Polygon[Fraction],
                       delta: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    for triangle in partition.triangular_map:
        neighbors = partition.triangle_neighbors(triangle)
        for other in partition.triangular_map:
            if other in neighbors:
                continue
            with pytest.raises(ValueError, match="don't share an edge"):
                partition.edge_between_triangles(triangle, other)


def test_sharing_only_a_vertex() -> None:
    partition, chunks = to_inner_triangle_partition()
    (bottom,) = chunks["bottom"].triangles
    (right,) = chunks["right"].triangles
    (inner,) = chunks["inner"].triangles

    assert partition.edge_between_triangles(bottom, inner) == Segment(
        to_point(1, 1), to_point(3, 1))
    # `bottom` and `right` share only the vertex (3, 1).
    with pytest.raises(ValueError, match="don't share an edge"):
        partition.edge_between_triangles(bottom, right)
