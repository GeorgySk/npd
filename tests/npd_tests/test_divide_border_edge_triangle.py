import math
from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import (border_edge_triangle,
                           divide_border_edge_triangle,
                           edge_neighbor)
from npd.structures import to_triangular_map
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   two_chunks_choices,
                                   varied_polygons)
from tests.utils import (border_length,
                         rebuilt_neighbors,
                         to_inner_triangle_partition,
                         to_point,
                         to_two_chunks,
                         unite_all)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       choices=two_chunks_choices)
def test_leaves_half_movable(polygon: Polygon[Fraction],
                             delta: Fraction,
                             operations: list[tuple[bool, int, int]],
                             choices: tuple[int, bool]) -> None:
    choice, in_remainder = choices
    setup = to_two_chunks(polygon, delta, operations, choice)
    if setup is None:
        return
    partition, chunk, neighbor = setup
    if in_remainder:
        chunk, neighbor = neighbor, chunk
    try:
        triangle = border_edge_triangle(chunk, neighbor, partition)
    except ValueError:
        return
    cut_neighbor = edge_neighbor(triangle, chunk, partition)
    (edge,) = (edge for edge, other in partition.triangular_map[
        triangle].items() if other is cut_neighbor)
    before = set(partition.triangular_map)
    areas = {chunk_: chunk_.area for chunk_ in partition.chunk_map}

    divide_border_edge_triangle(chunk, neighbor, partition)

    new = set(partition.triangular_map) - before
    assert before - set(partition.triangular_map) == {triangle,
                                                      cut_neighbor}
    assert all(edge.centroid in other.border.vertices for other in new)
    halves = {other for other in new
              if set(other.border.vertices)
              <= {*triangle.border.vertices, edge.centroid}}
    assert sum(half.area for half in halves) == triangle.area
    assert halves & set(partition.movable_triangles(chunk))
    assert partition.triangular_map == to_triangular_map(
        list(partition.triangular_map))
    assert partition.chunk_map == rebuilt_neighbors(partition)
    for chunk_ in partition.chunk_map:
        assert chunk_.area == areas[chunk_]
        assert chunk_.area == sum(other.area for other in chunk_.triangles)
        assert math.isclose(chunk_.perimeter,
                            border_length(chunk_.triangles, partition))


def test_inner_triangle() -> None:
    # Only `inner` is outside the chunk. `bottom` is the first triangle
    # with one neighbor in it and two in the chunk; it is cut at the
    # midpoint of its edge with `bottom_left`, the first of those two
    # that touches nothing else outside (`bottom_right` touches `right`,
    # which touches `inner`).
    partition, chunks = to_inner_triangle_partition()
    inner = chunks.pop("inner")
    (bottom,) = chunks["bottom"].triangles
    (bottom_left,) = chunks["bottom_left"].triangles
    chunk = unite_all(partition, *chunks.values())
    assert border_edge_triangle(chunk, inner, partition) is bottom

    divide_border_edge_triangle(chunk, inner, partition)

    assert bottom not in partition.triangular_map
    assert bottom_left not in partition.triangular_map
    midpoint = to_point(Fraction(5, 2), Fraction(1, 2))
    halves = [triangle for triangle in chunk.triangles
              if midpoint in triangle.border.vertices]
    assert len(halves) == 4  # noqa: PLR2004
    (inner_triangle,) = inner.triangles
    assert any(half in partition.movable_triangles(chunk)
               and inner_triangle in partition.triangle_neighbors(half)
               for half in halves)


def test_no_border_edge_triangle() -> None:
    # `inner` touches the chunk with all three edges, so no triangle of
    # the chunk has exactly one neighbor in it and two in the chunk.
    partition, chunks = to_inner_triangle_partition()
    inner = chunks.pop("inner")
    chunk = unite_all(partition, *chunks.values())

    with pytest.raises(ValueError, match="Couldn't find"):
        divide_border_edge_triangle(inner, chunk, partition)
