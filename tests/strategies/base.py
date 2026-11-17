from fractions import Fraction
from math import floor

from hypothesis import strategies as st
from hypothesis.strategies import SearchStrategy
from gon.base import (Contour,
                      Polygon,
                      Triangulation)
from ground.base import get_context
from hypothesis_geometry import planar

from npd.npd import steiner_points
from npd.structures import (CachedPolygon,
                            to_triangular_map)

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
# Generating Steiner points is cheap (no triangulation), so its tests
# can afford non-convex polygons with holes. Holes need a wider grid than
# `coordinates` gives (on it no holes are generated at all), and even then
# only ~6% of polygons get one unless at least one is required.
wide_coordinates = st.fractions(0, 4, max_denominator=2)
varied_polygons = st.one_of(
    planar.polygons(wide_coordinates,
                    max_size=8),
    planar.polygons(wide_coordinates,
                    max_size=8,
                    min_holes_size=1,
                    max_holes_size=2,
                    max_hole_size=4))
steiner_points_counts = st.integers(MIN_STEINER_POINTS_COUNT,
                                    MAX_STEINER_POINTS_COUNT)
unit_area_fractions = st.fractions(MIN_REQUIREMENT_FRACTION,
                                   MAX_REQUIREMENT_FRACTION)
deltas = st.fractions(Fraction(1, 10), Fraction(1))
# Steiner points spaced at least this far apart keep partitions to a few
# dozen triangles, so uniting them stays fast.
sparse_deltas = st.fractions(Fraction(1, 2), Fraction(1))
MAX_POLYGONS_LIST_SIZE = 4


def to_polygon_with_triangles(polygon: Polygon[Fraction],
                              delta: Fraction
                              ) -> tuple[Polygon[Fraction],
                                         list[CachedPolygon]]:
    """Triangulates the polygon the same way `Partition` does."""
    triangulation = Triangulation.constrained_delaunay(
        polygon,
        extra_points=steiner_points(polygon, delta),
        context=get_context())
    return polygon, list(map(CachedPolygon, triangulation.triangles()))


polygons_with_triangles = st.builds(to_polygon_with_triangles,
                                    varied_polygons,
                                    sparse_deltas)
triangles = planar.triangular_contours(wide_coordinates).map(CachedPolygon)
triangles_lists = st.lists(triangles, min_size=1, max_size=4)
scales = st.floats(1e-3, 1e3)
positive_fractions = st.fractions(0).filter(bool)
# Operations on a partition: (is a move rather than a unite, source chunk
# index for a unite, index of the chosen neighbor or move). Indices are
# taken modulo the available options.
MAX_PARTITION_OPERATIONS_COUNT = 10
partition_operations = st.lists(st.tuples(st.booleans(),
                                          st.integers(0),
                                          st.integers(0)),
                                max_size=MAX_PARTITION_OPERATIONS_COUNT)
# requirements as fractions of the area that `nonconvex.split` rejects
invalid_requirement_factors = (st.fractions(max_value=0)
                               | st.fractions(min_value=1))
triangular_maps = polygons_with_triangles.map(
    lambda polygon_with_triangles: to_triangular_map(
        polygon_with_triangles[1]))
triangles_lists_with_permutations = polygons_with_triangles.flatmap(
    lambda polygon_with_triangles: st.tuples(
        st.just(polygon_with_triangles[1]),
        st.permutations(polygon_with_triangles[1])))
# Fine-grained coordinates with signs, so lengths need rounding.
length_coordinates = st.fractions(-100, 100, max_denominator=1000)
length_segments = planar.segments(length_coordinates)
length_contours = planar.contours(length_coordinates)


def to_reoriented(contour: Contour[Fraction],
                  start: int,
                  reverse: bool) -> Contour[Fraction]:
    """Same contour starting at another vertex, maybe traversed backwards."""
    vertices = contour.vertices[start:] + contour.vertices[:start]
    return Contour(vertices[::-1] if reverse else vertices)


contours_with_reoriented = length_contours.flatmap(
    lambda contour: st.tuples(
        st.just(contour),
        st.builds(to_reoriented,
                  st.just(contour),
                  st.integers(0, len(contour.vertices) - 1),
                  st.booleans())))
multisegments = planar.multisegments(length_coordinates)
polygons_lists = st.lists(polygons,
                          min_size=1,
                          max_size=MAX_POLYGONS_LIST_SIZE)
polygons_lists_with_permutations = polygons_lists.flatmap(
    lambda polygons_list: st.tuples(st.just(polygons_list),
                                    st.permutations(polygons_list)))

MIN_PARTITION_SIZE = 1


def to_partitions(sum_: Fraction,  # noqa: C901
                  *,
                  min_value: Fraction = Fraction(0),
                  size: int = MIN_PARTITION_SIZE,
                  base: SearchStrategy[Fraction] | None = None
                  ) -> SearchStrategy[list[Fraction]]:
    if base is None:
        base = st.integers().map(Fraction)
    if size < MIN_PARTITION_SIZE:
        raise ValueError("`size` should not be less "
                         f"than {MIN_PARTITION_SIZE}.")
    if not (0 <= min_value <= sum_):
        raise ValueError(f"`min_value` should be in [0, {sum_}] interval.")
    if min_value:
        max_size = floor(sum_ / min_value)
        if max_size < size:
            raise ValueError(f"`size` should not be greater than {max_size}.")

    def to_proportions(numbers: list[Fraction]) -> list[Fraction]:
        return [2 * abs(number) / (1 + number * number) for number in numbers]

    def scale_to_sum(proportions: list[Fraction]) -> list[Fraction]:
        factor = sum_ / sum(proportions)
        return [proportion * factor for proportion in proportions]

    def bound_minimum(partition: list[Fraction]) -> list[Fraction]:
        minimum = min(partition)
        if minimum >= min_value:
            return partition
        partition_size = len(partition)
        denominator = sum_ - partition_size * minimum
        slope = sum_ - partition_size * min_value
        intercept = sum_ * (min_value - minimum)
        return [max((part * slope + intercept) / denominator, min_value)
                for part in partition]

    def normalize(partition: list[Fraction]) -> list[Fraction]:
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

    def is_valid(partition: list[Fraction]) -> bool:
        return sum(partition) == sum_

    return (st.lists(base,
                     min_size=size,
                     max_size=size)
            .filter(any)
            .map(to_proportions)
            .map(scale_to_sum)
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
