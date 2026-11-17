from npd.structures import Partition


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
