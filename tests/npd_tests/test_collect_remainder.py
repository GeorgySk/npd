import math
from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import (given,
                        strategies as st)

from npd.nonconvex import collect_remainder
from npd.npd import to_partition
from tests.strategies.base import (partition_operations,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (apply,
                         border_length,
                         rebuilt_neighbors,
                         to_strip_partition)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       choice=st.integers(0))
def test_unites_all_other_chunks(polygon: Polygon[Fraction],
                                 delta: Fraction,
                                 operations: list[tuple[bool, int, int]],
                                 choice: int) -> None:
    partition = to_partition(polygon, delta=delta)
    apply(partition, operations)
    if len(partition.chunk_map) < 2:  # noqa: PLR2004
        return
    candidates = [chunk for chunk in partition.chunk_map
                  if not partition.is_articulation_chunk(chunk)]
    chunk = candidates[choice % len(candidates)]
    triangles = dict(chunk.triangles)
    others = {triangle
              for other in partition.chunk_map
              if other is not chunk
              for triangle in other.triangles}

    result = collect_remainder(chunk, partition)

    assert len(partition.chunk_map) == 2  # noqa: PLR2004
    assert set(partition.chunk_map) == {chunk, result}
    assert chunk.triangles == triangles
    assert set(result.triangles) == others
    assert len(result.triangles) == len(others)
    assert result.area == partition.total_area - chunk.area
    assert math.isclose(result.perimeter,
                        border_length(result.triangles, partition))
    assert partition.chunk_map == rebuilt_neighbors(partition)


def test_rejects_articulation_chunk() -> None:
    partition, (_, middle, _) = to_strip_partition()

    with pytest.raises(ValueError, match="connected"):
        collect_remainder(middle, partition)
