import math
from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd.nonconvex import readjust
from npd.structures import to_triangular_map
from tests.strategies.base import (partition_operations,
                                   refinements,
                                   sparse_deltas,
                                   varied_polygons)
from tests.utils import (border_length,
                         is_connected,
                         rebuilt_neighbors,
                         to_two_chunks)


@given(polygon=varied_polygons,
       delta=sparse_deltas,
       operations=partition_operations,
       refinement=refinements)
def test_meets_requirement(polygon: Polygon[Fraction],
                           delta: Fraction,
                           operations: list[tuple[bool, int, int]],
                           refinement: tuple[int, Fraction]) -> None:
    choice, fraction = refinement
    setup = to_two_chunks(polygon, delta, operations, choice)
    if setup is None:
        return
    partition, chunk, remainder = setup
    requirement = fraction * chunk.area
    connected = [chunk_ for chunk_ in (chunk, remainder)
                 if is_connected(partition, chunk_.triangles)]

    result = readjust(chunk,
                      remainder,
                      requirement=requirement,
                      partition=partition)

    assert result == (chunk, remainder)
    assert chunk.area == requirement
    assert remainder.area == partition.total_area - requirement
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
def test_rejects_requirement_not_below_area(polygon: Polygon[Fraction],
                                            delta: Fraction,
                                            operations: list[tuple[bool,
                                                                   int,
                                                                   int]],
                                            refinement: tuple[int, Fraction]
                                            ) -> None:
    choice, fraction = refinement
    setup = to_two_chunks(polygon, delta, operations, choice)
    if setup is None:
        return
    partition, chunk, remainder = setup

    for requirement in (chunk.area, chunk.area * (1 + fraction)):
        with pytest.raises(ValueError, match="below the chunk's area"):
            readjust(chunk,
                     remainder,
                     requirement=requirement,
                     partition=partition)
