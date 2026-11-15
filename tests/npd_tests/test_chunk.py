import math
from fractions import Fraction

from hypothesis import given

from npd.structures import (CachedPolygon,
                            Chunk,
                            to_chunk)
from tests.strategies.base import (scales,
                                   triangles,
                                   triangles_lists)

# The isoperimetric inequality: a circle has the largest area for its
# perimeter, and its compactness is sqrt(pi r^2) / (2 pi r).
MAX_COMPACTNESS = 1 / (2 * math.sqrt(math.pi))


@given(triangles_lists)
def test_hash_survives_triangle_changes(
        triangles_list: list[CachedPolygon]) -> None:
    first, *rest = triangles_list
    chunk = to_chunk(first)
    chunks = {chunk: None}

    for triangle in rest:
        chunk.triangles.setdefault(triangle)
    chunk.triangles.pop(first)

    assert chunk in chunks


@given(triangles)
def test_equality_is_identity(triangle: CachedPolygon) -> None:
    assert to_chunk(triangle) != to_chunk(triangle)


@given(triangles)
def test_compactness_bounds(triangle: CachedPolygon) -> None:
    assert 0 < to_chunk(triangle).compactness() <= MAX_COMPACTNESS


@given(triangles, scales)
def test_compactness_scale_invariance(triangle: CachedPolygon,
                                      scale: float) -> None:
    chunk = to_chunk(triangle)
    scaled = Chunk(triangles=chunk.triangles,
                   area=chunk.area * Fraction(scale) ** 2,
                   perimeter=chunk.perimeter * scale)

    assert math.isclose(scaled.compactness(), chunk.compactness())
