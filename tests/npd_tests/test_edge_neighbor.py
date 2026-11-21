from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import edge_neighbor
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         to_inner_triangle_partition,
                         to_point,
                         to_polygon,
                         unite_all)


def touches_other_chunk(triangle: CachedPolygon,
                        chunk: Chunk,
                        partition: Partition) -> bool:
    return any(other is not None and other not in chunk.triangles
               for other in partition.triangular_map[triangle].values())


def meets_condition(candidate: CachedPolygon,
                    triangle: CachedPolygon,
                    chunk: Chunk,
                    partition: Partition) -> bool:
    """
    The condition as stated in the docstring.

    The candidate reaches the chunk's border only through the triangle:
    none of its other neighbors is outside the chunk or touches it.
    """
    return not any(
        other not in chunk.triangles
        or touches_other_chunk(other, chunk, partition)
        for other in partition.triangular_map[candidate].values()
        if other is not None and other is not triangle)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_prefers_first_qualifying_neighbor(polygon: Polygon[Fraction],
                                           delta: Fraction,
                                           operations: list[tuple[bool,
                                                                  int,
                                                                  int]]
                                           ) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        for triangle in chunk.triangles:
            if not touches_other_chunk(triangle, chunk, partition):
                continue
            candidates = [candidate
                          for candidate in partition.triangle_neighbors(
                              triangle)
                          if candidate in chunk.triangles]
            if not candidates:
                with pytest.raises(ValueError, match="Couldn't find"):
                    edge_neighbor(triangle, chunk, partition)
                continue
            expected = next(
                (candidate for candidate in candidates
                 if meets_condition(candidate, triangle, chunk, partition)),
                candidates[0])
            assert edge_neighbor(triangle, chunk, partition) is expected


def test_skips_neighbor_near_border() -> None:
    partition, chunks = to_inner_triangle_partition()
    (bottom,) = chunks["bottom"].triangles
    (bottom_left,) = chunks["bottom_left"].triangles
    chunk = unite_all(partition, *(chunk for name, chunk in chunks.items()
                                   if name != "bottom_right"))

    # `inner` is skipped: its neighbor `right` touches `bottom_right`.
    assert edge_neighbor(bottom, chunk, partition) is bottom_left


def test_falls_back_to_first_chunk_neighbor() -> None:
    partition, chunks = to_inner_triangle_partition()
    (inner,) = chunks["inner"].triangles
    chunk = unite_all(partition,
                      chunks["inner"], chunks["bottom"], chunks["left"])
    # both `bottom` and `left` touch other chunks besides through `inner`
    first = next(neighbor for neighbor in partition.triangle_neighbors(inner)
                 if neighbor in chunk.triangles)

    assert edge_neighbor(inner, chunk, partition) is first


def test_skips_neighbor_touching_other_chunk_directly() -> None:
    # Only `outside` (around the inner vertex (2, 1) with `triangle` and
    # `first`) is in the other chunk. `first` touches it directly, but
    # as all of `outside`'s neighbors are in the chunk, `first` used to
    # be taken as reaching the chunk's border only through `triangle`.
    partition = Partition(to_polygon([(0, 0), (6, 0), (7, 5), (3, 5)]),
                          extra_points=[to_point(2, 1)])
    chunks = {frozenset((int(vertex.x), int(vertex.y))
                        for vertex in chunk_triangle.border.vertices): chunk
              for chunk in partition.chunk_map
              for chunk_triangle in chunk.triangles}
    (triangle,) = chunks[frozenset({(2, 1), (6, 0), (3, 5)})].triangles
    (second,) = chunks[frozenset({(3, 5), (6, 0), (7, 5)})].triangles
    chunk = unite_all(partition, *(
        chunk for vertices, chunk in chunks.items()
        if vertices != frozenset({(0, 0), (2, 1), (3, 5)})))
    first = next(neighbor for neighbor in partition.triangle_neighbors(
        triangle) if neighbor in chunk.triangles)
    assert first is not second

    assert edge_neighbor(triangle, chunk, partition) is second


def test_no_chunk_neighbor() -> None:
    partition, chunks = to_inner_triangle_partition()
    (inner,) = chunks["inner"].triangles

    with pytest.raises(ValueError, match="Couldn't find"):
        edge_neighbor(inner, chunks["inner"], partition)
