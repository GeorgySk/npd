import math
from fractions import Fraction
from math import floor
from numbers import Real
from typing import List

from hypothesis import strategies as st
from hypothesis.strategies import SearchStrategy
from hypothesis_geometry import planar

from npd.npd import (Partition,
                     to_partition)
from tests.config import (MAX_COORDINATE,
                          MAX_DENOMINATOR,
                          MAX_HOLES_SIZE,
                          MAX_HOLE_SIZE,
                          MAX_POLYGON_SIZE,
                          MAX_REQRUIREMENTS_COUNT,
                          MAX_REQUIREMENT_FRACTION,
                          MAX_STEINER_POINTS_COUNT,
                          MIN_COORDINATE,
                          MIN_HOLES_SIZE,
                          MIN_HOLE_SIZE,
                          MIN_POLYGON_SIZE,
                          MIN_REQRUIREMENTS_COUNT,
                          MIN_REQUIREMENT_FRACTION,
                          MIN_STEINER_POINTS_COUNT)

coordinates = st.fractions(MIN_COORDINATE, MAX_COORDINATE,
                           max_denominator=MAX_DENOMINATOR)
polygons = planar.polygons(coordinates,
                           min_size=MIN_POLYGON_SIZE,
                           max_size=MAX_POLYGON_SIZE,
                           min_holes_size=MIN_HOLES_SIZE,
                           max_holes_size=MAX_HOLES_SIZE,
                           min_hole_size=MIN_HOLE_SIZE,
                           max_hole_size=MAX_HOLE_SIZE)
steiner_points_counts = st.integers(MIN_STEINER_POINTS_COUNT,
                                    MAX_STEINER_POINTS_COUNT)
unit_area_fractions = st.fractions(MIN_REQUIREMENT_FRACTION,
                                   MAX_REQUIREMENT_FRACTION)
partitions = st.builds(to_partition,
                       polygon=polygons,
                       steiner_points_count=steiner_points_counts)

MIN_PARTITION_SIZE = 1


def to_partitions(sum_: Real,
                  *,
                  min_value: Real = 0,
                  size: int = MIN_PARTITION_SIZE,
                  base: SearchStrategy[Real] = st.integers()
                  ) -> SearchStrategy[List[Real]]:
    if size < MIN_PARTITION_SIZE:
        raise ValueError('`size` should not be less '
                         f'than {MIN_PARTITION_SIZE}.')
    if not (0 <= min_value <= sum_):
        raise ValueError(f'`min_value` should be in [0, {sum_}] interval.')
    if min_value:
        max_size = floor(sum_ / min_value)
        if max_size < size:
            raise ValueError(f'`size` should not be greater than {max_size}.')

    def to_proportions(numbers: List[Real]) -> List[Real]:
        return [2 * abs(number) / (1 + number * number) for number in numbers]

    def to_partition(proportions: List[Real]) -> List[Real]:
        factor = sum_ / sum(proportions)
        return [proportion * factor for proportion in proportions]

    def bound_minimum(partition: List[Real]) -> List[Real]:
        minimum = min(partition)
        if minimum >= min_value:
            return partition
        partition_size = len(partition)
        denominator = sum_ - partition_size * minimum
        slope = sum_ - partition_size * min_value
        intercept = sum_ * (min_value - minimum)
        return [max((part * slope + intercept) / denominator, min_value)
                for part in partition]

    def normalize(partition: List[Real]) -> List[Real]:
        partition_sum = sum(partition)
        if partition_sum < sum_:
            arg_min = min(range(len(partition)),
                          key=partition.__getitem__)
            partition[arg_min] += sum_ - partition_sum
        elif partition_sum > sum_:
            arg_max = max(range(len(partition)),
                          key=partition.__getitem__)
            partition[arg_max] -= partition_sum - sum_
        return partition

    def is_valid(partition: List[Real]) -> bool:
        return sum(partition) == sum_

    return (st.lists(base,
                     min_size=size,
                     max_size=size)
            .filter(any)
            .map(to_proportions)
            .map(to_partition)
            .map(bound_minimum)
            .map(normalize)
            .filter(is_valid))


requirements_fractions = st.integers(MIN_REQRUIREMENTS_COUNT,
                                     MAX_REQRUIREMENTS_COUNT
                                     ).flatmap(
    lambda n: to_partitions(Fraction(1),
                            base=st.fractions(),
                            min_value=MIN_REQUIREMENT_FRACTION,
                            size=n))
