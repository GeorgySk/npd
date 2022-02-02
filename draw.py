import inspect
import random
from functools import singledispatch
from typing import (List,
                    Union)

import matplotlib.pyplot as plt
from gon.base import (Contour,
                      Geometry,
                      Point,
                      Polygon,
                      Segment)

from npd.structures import Chunk


def draw(*geometries: Union[Geometry, Chunk, List], **kwargs) -> None:
    if not draw.DRAW:
        return
    colored_geometries = {color: kwargs.pop(color, [])
                          for color in ['red', 'blue']}
    for geometry in geometries:
        _draw(geometry, color='black', **kwargs)
    for color, geometries in colored_geometries.items():
        if not geometries:
            continue
        if not isinstance(geometries, list):
            geometries = [geometries]
        for geometry in geometries:
            _draw(geometry, color=color)
    caller = inspect.stack()[1]
    plt.title(f"{caller.function}:{caller.lineno}")
    plt.axis('scaled')
    plt.show(block=True)


@singledispatch
def _draw(geometry: Union[Geometry, Chunk, list], *args, **kwargs
          ) -> None:
    raise TypeError(f"Got unexpected geometry type: {type(geometry)}")


@_draw.register(Point)
def _(geometry: Point, *args, **kwargs) -> None:
    kwargs.setdefault('marker', 'o')
    plt.plot(geometry.x, geometry.y, *args, **kwargs)


@_draw.register(Segment)
def _(segment: Segment, *args, label: str = None, **kwargs):
    start, end = segment.start, segment.end
    _draw(start, **({'color': kwargs['color']}
                    if 'color' in kwargs else {}))
    plt.plot(*zip((start.x, start.y), (end.x, end.y)), *args, **kwargs)


@_draw.register(Contour)
def _(geometry: Contour, *args, fill: bool = False, **kwargs) -> None:
    vertices = [(point.x, point.y) for point in geometry.vertices]
    xs, ys = zip(*(vertices + vertices[:1]))
    xs = list(map(float, xs))
    ys = list(map(float, ys))
    if fill:
        plt.fill(xs, ys, *args, **kwargs)
    plt.plot(xs, ys, *args, **kwargs)


@_draw.register(Polygon)
def _(geometry: Polygon, *args, **kwargs) -> None:
    border_color = kwargs.pop('color')
    fill_color = kwargs.pop('fill_color', None)
    if fill_color is not None:
        _draw(geometry.border, *args, fill=True, color=fill_color, **kwargs)
    _draw(geometry.border, *args, color=border_color, **kwargs)
    for hole in geometry.holes:
        _draw(hole, *args, fill=True, color=(0, 0, 0, 0.1), **kwargs)


@_draw.register(Chunk)
def _(geometry: Chunk, *args, **kwargs) -> None:
    random_color = (random.uniform(0, 1),
                    random.uniform(0, 1),
                    random.uniform(0, 1),
                    0.8)
    for triangle in geometry.triangles:
        _draw(triangle, *args, fill_color=random_color, **kwargs)


@_draw.register(list)
def _(geometry: list, *args, **kwargs) -> None:
    for i, point in enumerate(geometry):
        _draw(point, *args, **kwargs)
        plt.text(point.x + 0.05, point.y + 0.05, str(i))
