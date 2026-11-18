from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import border_edge_triangle
from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         to_inner_triangle_partition,
                         unite_all)


def meets_conditions(triangle: CachedPolygon,
                     chunk: Chunk,
                     neighbor: Chunk,
                     partition: Partition) -> bool:
    """The conditions as stated in the docstring, worked out per edge."""
    owners = {other: owner
              for owner in partition.chunk_map
              for other in owner.triangles}
    across = list(partition.triangular_map[triangle].values())
    in_neighbor = [other for other in across
                   if other is not None and owners[other] is neighbor]
    in_chunk = [other for other in across
                if other is not None and owners[other] is chunk]
    if len(in_neighbor) != 1 or len(in_chunk) != 2:  # noqa: PLR2004
        return False
    vertices = set(triangle.border.vertices)
    touching = [
        other for other in neighbor.triangles
        if other is not in_neighbor[0]
        and vertices & set(other.border.vertices)
        and any(beyond is not None and owners[beyond] is not neighbor
                for beyond in partition.triangular_map[other].values())]
    return len(touching) <= 1


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_returns_first_qualifying_triangle(polygon: Polygon[Fraction],
                                           delta: Fraction,
                                           operations: list[tuple[bool,
                                                                  int,
                                                                  int]]
                                           ) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        for neighbor in partition.chunk_neighbors(chunk):
            expected = next(
                (triangle for triangle in chunk.triangles
                 if meets_conditions(triangle, chunk, neighbor, partition)),
                None)
            if expected is None:
                with pytest.raises(ValueError, match="Couldn't find"):
                    border_edge_triangle(chunk, neighbor, partition)
            else:
                assert border_edge_triangle(chunk,
                                            neighbor,
                                            partition) is expected


@pytest.mark.parametrize(("neighbor_names", "found"),
                         [(["right"], True),
                          # one triangle besides `right` sharing a vertex
                          # with `inner` and touching the chunk
                          (["right", "bottom_right"], True),
                          # two such triangles
                          (["right", "bottom_right", "top"], False)])
def test_inner_triangle(neighbor_names: list[str], found: bool) -> None:
    partition, chunks = to_inner_triangle_partition()
    (inner,) = chunks["inner"].triangles
    chunk = unite_all(partition,
                      chunks["inner"], chunks["bottom"], chunks["left"])
    neighbor = unite_all(partition,
                         *(chunks[name] for name in neighbor_names))

    if found:
        assert border_edge_triangle(chunk, neighbor, partition) is inner
    else:
        with pytest.raises(ValueError, match="Couldn't find"):
            border_edge_triangle(chunk, neighbor, partition)


def test_no_triangle_with_one_neighbor_in_neighbor_chunk() -> None:
    partition, chunks = to_inner_triangle_partition()
    chunk = unite_all(partition,
                      chunks["inner"], chunks["bottom"], chunks["left"])

    # `bottom` has one neighbor in it but only one in its own chunk
    with pytest.raises(ValueError, match="Couldn't find"):
        border_edge_triangle(chunk, chunks["bottom_right"], partition)
