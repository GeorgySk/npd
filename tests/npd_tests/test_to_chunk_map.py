import math

from hypothesis import given

from npd.structures import (Chunk,
                            CachedPolygon,
                            TriangularMap,
                            to_chunk_map)
from tests.strategies.base import triangular_maps


def to_triangle(chunk: Chunk) -> CachedPolygon:
    (triangle,) = chunk.triangles
    return triangle


@given(triangular_maps)
def test_one_chunk_per_triangle(triangular_map: TriangularMap) -> None:
    result = to_chunk_map(triangular_map)

    assert all(len(chunk.triangles) == 1 for chunk in result)
    assert list(map(to_triangle, result)) == list(triangular_map)


@given(triangular_maps)
def test_chunk_measures(triangular_map: TriangularMap) -> None:
    result = to_chunk_map(triangular_map)

    for chunk in result:
        triangle = to_triangle(chunk)
        assert chunk.area == triangle.area
        assert math.isclose(chunk.perimeter,
                             sum(math.dist((edge.start.x, edge.start.y),
                                           (edge.end.x, edge.end.y))
                                 for edge in triangle.edges))


@given(triangular_maps)
def test_neighbors_match_triangular_map(triangular_map: TriangularMap
                                        ) -> None:
    result = to_chunk_map(triangular_map)

    for chunk, neighbors in result.items():
        expected = {frozenset([edge]): neighbor
                    for edge, neighbor
                    in triangular_map[to_triangle(chunk)].items()
                    if neighbor is not None}
        assert {edges: to_triangle(neighbor)
                for edges, neighbor in neighbors.items()} == expected


@given(triangular_maps)
def test_neighbors_are_keys(triangular_map: TriangularMap) -> None:
    result = to_chunk_map(triangular_map)

    assert all(neighbor in result and result[neighbor][edges] is chunk
               for chunk, neighbors in result.items()
               for edges, neighbor in neighbors.items())

