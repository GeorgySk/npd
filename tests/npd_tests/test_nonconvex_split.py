from fractions import Fraction

import pytest
from gon.base import Polygon
from hypothesis import given

from npd import nonconvex
from npd.npd import to_partition
from tests.strategies.base import (invalid_requirement_factors,
                                   polygons,
                                   sparse_deltas)


@given(polygon=polygons,
       delta=sparse_deltas,
       factor=invalid_requirement_factors)
def test_invalid_requirement(polygon: Polygon[Fraction],
                             delta: Fraction,
                             factor: Fraction) -> None:
    partition = to_partition(polygon, delta=delta)

    with pytest.raises(ValueError, match="positive and less than"):
        nonconvex.split(partition, requirement=polygon.area * factor)
