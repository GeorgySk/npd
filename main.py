import logging
from fractions import Fraction

from gon.base import (Contour,
                      Point,
                      Polygon)
from ground.base import get_context

from draw import draw
from npd.npd import split

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("matplotlib.font_manager").disabled = True
context = get_context()

draw.DRAW = True


if __name__ == '__main__':
    polygon = Polygon(Contour([Point(Fraction(0, 1), Fraction(0, 1)),
                               Point(Fraction(2, 1), Fraction(0, 1)),
                               Point(Fraction(1, 2), Fraction(1, 1)),
                               Point(Fraction(1, 2), Fraction(2, 1))]), [])
    unit_area_fraction = Fraction(1, 2)
    steiner_points_count = 3
    requirements_fractions = [unit_area_fraction, 1 - unit_area_fraction]
    requirements = [requirement * polygon.area
                    for requirement in requirements_fractions]
    parts = split(polygon,
                  requirements=requirements,
                  steiner_points_count=steiner_points_count)
    draw(polygon, *parts)
