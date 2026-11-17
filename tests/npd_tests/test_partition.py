import math
from collections.abc import Collection
from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon,
                      Segment)
from hypothesis import given
import pytest

from npd.npd import to_partition
from npd.structures import (CachedPolygon,
                            Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)


def apply(partition: Partition,
          operations: list[tuple[bool, int, int]]) -> None:
    for is_move, source_index, choice_index in operations:
        if is_move:
            owners = {triangle: chunk
                      for chunk in partition.chunk_map
                      for triangle in chunk.triangles}
            # a move must leave its source non-empty
            moves = [(triangle, source, owners[neighbor])
                     for source in partition.chunk_map
                     if len(source.triangles) > 1
                     for triangle in partition.movable_triangles(source)
                     for neighbor in partition.triangle_neighbors(triangle)
                     if owners[neighbor] is not source]
            if not moves:
                continue
            triangle, source, target = moves[choice_index % len(moves)]
            partition.move(triangle, source=source, target=target)
        else:
            chunks = list(partition.chunk_map)
            source = chunks[source_index % len(chunks)]
            neighbors = partition.chunk_neighbors(source)
            if not neighbors:
                continue
            partition.unite(source,
                            neighbors[choice_index % len(neighbors)])


def rebuilt_neighbors(partition: Partition
                      ) -> dict[Chunk, dict[Chunk,
                                            frozenset[Segment[Fraction]]]]:
    owners: dict[CachedPolygon, Chunk] = {
        triangle: chunk
        for chunk in partition.chunk_map
        for triangle in chunk.triangles}
    result: dict[Chunk, dict[Chunk, set[Segment[Fraction]]]] = {
        chunk: {} for chunk in partition.chunk_map}
    for triangle, neighbors in partition.triangular_map.items():
        chunk = owners[triangle]
        for edge, neighbor in neighbors.items():
            if neighbor is not None and owners[neighbor] is not chunk:
                result[chunk].setdefault(owners[neighbor], set()).add(edge)
    return {chunk: {neighbor: frozenset(edges)
                    for neighbor, edges in neighbors_edges.items()}
            for chunk, neighbors_edges in result.items()}


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_chunk_map_consistency(polygon: Polygon[Fraction],
                               delta: Fraction,
                               operations: list[tuple[bool, int, int]]
                               ) -> None:
    partition = to_partition(polygon, delta=delta)

    apply(partition, operations)

    assert partition.chunk_map == rebuilt_neighbors(partition)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_chunks_partition_triangles(polygon: Polygon[Fraction],
                                    delta: Fraction,
                                    operations: list[tuple[bool, int, int]]
                                    ) -> None:
    partition = to_partition(polygon, delta=delta)

    apply(partition, operations)

    owned = [triangle
             for chunk in partition.chunk_map
             for triangle in chunk.triangles]
    assert all(chunk.triangles for chunk in partition.chunk_map)
    assert len(owned) == len(partition.triangular_map)
    assert set(owned) == set(partition.triangular_map)


def test_unite_merges_common_neighbor_edges() -> None:
    # Three triangles around the centroid, each touching the other two.
    border = Contour([Point(Fraction(0), Fraction(0)),
                      Point(Fraction(3), Fraction(0)),
                      Point(Fraction(0), Fraction(3))])
    partition = Partition(Polygon(border, []),
                          extra_points=[Point(Fraction(1), Fraction(1))])
    first, second, third = partition.chunk_map
    first_edges = partition.shared_edges(third, first)
    second_edges = partition.shared_edges(third, second)
    assert first_edges
    assert second_edges

    union = partition.unite(first, second)

    assert partition.chunk_map[third] == {union: first_edges | second_edges}
    assert partition.chunk_map == rebuilt_neighbors(partition)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_shared_edges(polygon: Polygon[Fraction],
                      delta: Fraction,
                      operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    expected = rebuilt_neighbors(partition)
    assert all(partition.shared_edges(chunk, other)
               == expected[chunk].get(other, frozenset())
               for chunk in partition.chunk_map
               for other in partition.chunk_map)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_unite_rejects_non_neighbors(polygon: Polygon[Fraction],
                                     delta: Fraction,
                                     operations: list[tuple[bool, int, int]]
                                     ) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    chunk_map = {chunk: dict(neighbors)
                 for chunk, neighbors in partition.chunk_map.items()}

    for chunk in chunk_map:
        for other in chunk_map:
            if other not in partition.chunk_neighbors(chunk):
                with pytest.raises(ValueError, match="neighbor"):
                    partition.unite(chunk, other)
    assert partition.chunk_map == chunk_map


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_chunk_measures(polygon: Polygon[Fraction],
                        delta: Fraction,
                        operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        assert chunk.area == sum((triangle.area
                                  for triangle in chunk.triangles),
                                 Fraction(0))
        assert math.isclose(chunk.perimeter,
                            border_length(chunk.triangles, partition))


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_union_perimeter(polygon: Polygon[Fraction],
                         delta: Fraction,
                         operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)

    for chunk in partition.chunk_map:
        for other in partition.chunk_map:
            if other is chunk:
                continue
            result = partition.union_perimeter(chunk, other)
            assert math.isclose(
                result,
                border_length({**chunk.triangles, **other.triangles},
                              partition))
            assert math.isclose(result,
                                partition.union_perimeter(other, chunk))


def border_length(triangles: Collection[CachedPolygon],
                  partition: Partition) -> float:
    return sum(math.dist((edge.start.x, edge.start.y),
                         (edge.end.x, edge.end.y))
               for triangle in triangles
               for edge, neighbor in partition.triangular_map[triangle].items()
               if neighbor not in triangles)
