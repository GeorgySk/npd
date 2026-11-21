import math
from fractions import Fraction

from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import most_compact_pair
from npd.npd import to_partition
from npd.structures import (Chunk,
                            Partition)
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         border_length,
                         to_strip_partition)


def scores(partition: Partition) -> dict[Chunk, float]:
    """Scores eligible chunks, measuring borders from the triangles."""
    result = {}
    for chunk in partition.chunk_map:
        if partition.is_articulation_chunk(chunk):
            continue
        rest = {triangle: None
                for other in partition.chunk_map
                if other is not chunk
                for triangle in other.triangles}
        chunk_compactness = (chunk.area ** 0.5
                             / border_length(chunk.triangles, partition))
        rest_compactness = ((partition.total_area - chunk.area) ** 0.5
                            / border_length(rest, partition))
        result[chunk] = min(chunk_compactness, rest_compactness)
    return result


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations)
def test_picks_best_chunk(polygon: Polygon[Fraction],
                          delta: Fraction,
                          operations: list[tuple[bool, int, int]]) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    if len(partition.chunk_map) < 2:  # noqa: PLR2004
        return
    expected = scores(partition)
    triangles = set(partition.triangular_map)

    chunk, remainder = most_compact_pair(partition)

    assert chunk in expected
    assert math.isclose(expected[chunk], max(expected.values()))
    assert set(partition.chunk_map) == {chunk, remainder}
    assert set(remainder.triangles) == triangles - set(chunk.triangles)


def test_skips_articulation_chunk() -> None:
    partition, (left, middle, right) = to_strip_partition()
    sides = {*left.triangles, *right.triangles}

    chunk, remainder = most_compact_pair(partition)

    assert chunk is not middle
    assert set(chunk.triangles) <= sides
    assert set(remainder.triangles) >= set(middle.triangles)
