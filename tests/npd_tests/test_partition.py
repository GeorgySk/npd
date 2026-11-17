from fractions import Fraction

from gon.base import (Polygon,
                      Segment)
from hypothesis import given

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

    stored = {chunk: {neighbor: edges for edges, neighbor in neighbors.items()}
              for chunk, neighbors in partition.chunk_map.items()}
    assert all(len(stored[chunk]) == len(neighbors)
               for chunk, neighbors in partition.chunk_map.items())
    assert stored == rebuilt_neighbors(partition)
