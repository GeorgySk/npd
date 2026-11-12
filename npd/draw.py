"""Plotting of intermediate steps for debugging the algorithm."""
import inspect
import random
from functools import singledispatch

import matplotlib.pyplot as plt
from gon.base import (Contour,
                      Geometry,
                      Point,
                      Polygon,
                      Segment)

from npd.structures import Chunk

Color = str | tuple[float, float, float, float]
Drawable = Geometry | Chunk | list


def draw(*geometries: Drawable,
         red: Drawable | None = None,
         blue: Drawable | None = None) -> None:
    """
    Plots the geometries and blocks until the plot window is closed.

    Does nothing unless `draw.DRAW` is set to True.
    :param geometries: geometries to plot in black
    :param red: geometries to highlight in red
    :param blue: geometries to highlight in blue
    """
    if not draw.DRAW:
        return
    for geometry in geometries:
        _draw(geometry, color="black")
    for color, highlighted in (("red", red), ("blue", blue)):
        if not highlighted:
            continue
        for geometry in (highlighted if isinstance(highlighted, list)
                         else [highlighted]):
            _draw(geometry, color=color)
    caller = inspect.stack()[1]
    plt.title(f"{caller.function}:{caller.lineno}")
    plt.axis("scaled")
    plt.show(block=True)


@singledispatch
def _draw(geometry: Drawable, *, color: Color) -> None:  # noqa: ARG001
    raise TypeError(f"Got unexpected geometry type: {type(geometry)}")


@_draw.register(Point)
def _(geometry: Point, *, color: Color) -> None:
    plt.plot(geometry.x, geometry.y, marker="o", color=color)


@_draw.register(Segment)
def _(segment: Segment, *, color: Color) -> None:
    start, end = segment.start, segment.end
    _draw(start, color=color)
    plt.plot(*zip((start.x, start.y), (end.x, end.y)), color=color)


@_draw.register(Contour)
def _(geometry: Contour, *, color: Color) -> None:
    _draw_contour(geometry, color=color)


@_draw.register(Polygon)
def _(geometry: Polygon, *, color: Color) -> None:
    _draw_polygon(geometry, color=color)


@_draw.register(Chunk)
def _(geometry: Chunk, *, color: Color) -> None:
    random_color = (random.uniform(0, 1),
                    random.uniform(0, 1),
                    random.uniform(0, 1),
                    0.8)
    for triangle in geometry.triangles:
        _draw_polygon(triangle, color=color, fill_color=random_color)


@_draw.register(list)
def _(geometry: list, *, color: Color) -> None:
    for point in geometry:
        _draw(point, color=color)


def _draw_contour(contour: Contour,
                  *,
                  color: Color,
                  fill: bool = False) -> None:
    vertices = [(point.x, point.y) for point in contour.vertices]
    xs, ys = zip(*(vertices + vertices[:1]))
    xs = list(map(float, xs))
    ys = list(map(float, ys))
    if fill:
        plt.fill(xs, ys, color=color)
    plt.plot(xs, ys, color=color)


def _draw_polygon(polygon: Polygon,
                  *,
                  color: Color,
                  fill_color: Color | None = None) -> None:
    if fill_color is not None:
        _draw_contour(polygon.border, fill=True, color=fill_color)
    _draw_contour(polygon.border, color=color)
    for hole in polygon.holes:
        _draw_contour(hole, fill=True, color=(0, 0, 0, 0.1))
