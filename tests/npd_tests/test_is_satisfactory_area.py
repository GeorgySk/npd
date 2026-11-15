from fractions import Fraction

from hypothesis import given

from npd.nonconvex import is_satisfactory_area
from tests.strategies.base import positive_fractions


@given(area=positive_fractions,
       requirement=positive_fractions)
def test_stopping_rule(area: Fraction,
                       requirement: Fraction) -> None:
    # Stop once one more doubling would overshoot and leave a gap at
    # least as large as the current one.
    expected = (abs(requirement - area) <= abs(2 * area - requirement)
                and 2 * area >= requirement)

    assert is_satisfactory_area(area, requirement=requirement) is expected


@given(requirement=positive_fractions,
       shortfall=positive_fractions)
def test_boundary(requirement: Fraction,
                  shortfall: Fraction) -> None:
    threshold = Fraction(2, 3) * requirement

    assert is_satisfactory_area(threshold, requirement=requirement)
    assert not is_satisfactory_area(threshold - shortfall,
                                    requirement=requirement)


@given(area=positive_fractions,
       requirement=positive_fractions,
       increment=positive_fractions)
def test_monotonic(area: Fraction,
                   requirement: Fraction,
                   increment: Fraction) -> None:
    if is_satisfactory_area(area, requirement=requirement):
        assert is_satisfactory_area(area + increment,
                                    requirement=requirement)


@given(area=positive_fractions,
       requirement=positive_fractions,
       scale=positive_fractions)
def test_scale_invariance(area: Fraction,
                          requirement: Fraction,
                          scale: Fraction) -> None:
    assert (is_satisfactory_area(area * scale,
                                 requirement=requirement * scale)
            is is_satisfactory_area(area, requirement=requirement))
