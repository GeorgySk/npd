import math
from fractions import Fraction

from gon.base import (Contour,
                      Multisegment,
                      Segment)
from hypothesis import given
import pytest

from npd.utils import fast_length
from tests.strategies.base import (contours_with_reoriented,
                                   length_contours,
                                   length_segments,
                                   multisegments)


@given(length_segments)
def test_segment(segment: Segment[Fraction]) -> None:
    start, end = segment.start, segment.end

    assert math.isclose(fast_length(segment),
                        math.dist((float(start.x), float(start.y)),
                                  (float(end.x), float(end.y))))


@given(length_segments)
def test_segment_direction(segment: Segment[Fraction]) -> None:
    assert (fast_length(Segment(segment.end, segment.start))
            == fast_length(segment))


@given(length_contours)
def test_contour(contour: Contour[Fraction]) -> None:
    assert math.isclose(fast_length(contour),
                        math.fsum(map(fast_length, contour.segments)))


@given(contours_with_reoriented)
def test_contour_orientation(contour_with_reoriented: tuple[Contour[Fraction],
                                                            Contour[Fraction]]
                             ) -> None:
    contour, reoriented = contour_with_reoriented

    assert math.isclose(fast_length(reoriented), fast_length(contour))


@given(multisegments)
def test_unsupported(multisegment: Multisegment[Fraction]) -> None:
    with pytest.raises(TypeError, match="Cannot calculate length"):
        fast_length(multisegment)
