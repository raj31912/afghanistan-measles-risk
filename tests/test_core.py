import numpy as np

from measles_analysis.core import (
    P,
    age_eligibility_proxy,
    coverage_adjusted_transmission,
    preferential_sia,
    random_sia,
    routine_mcv1_gap,
    stable_province_seed,
)


def test_central_threshold_and_coverage_requirements():
    assert np.isclose(P.hit, 1 - 1 / 15)
    assert np.isclose(P.hit / P.ve_one_dose * 100, 100.35842293906809)
    assert np.isclose(P.hit / P.ve_two_dose * 100, 96.21993127147766)


def test_rcov_formula():
    assert np.isclose(coverage_adjusted_transmission(0.512), 7.8576)


def test_preferential_and_random_sia_are_distinct():
    c = 0.512
    p_immune = c * P.ve_one_dose
    expected_pref = max((P.hit - p_immune) / P.ve_one_dose, 0)
    expected_random = max((P.hit - p_immune) / ((1 - p_immune) * P.ve_one_dose), 0)
    assert np.isclose(preferential_sia(c), expected_pref)
    assert np.isclose(random_sia(c), expected_random)
    assert random_sia(c) > preferential_sia(c)


def test_sia_is_not_clipped_at_100_percent():
    # At sufficiently low coverage, the theoretical random-targeting requirement
    # can exceed 1.0. The analysis intentionally retains this infeasibility signal.
    x = random_sia(0.0, r0=18, ve=0.93)
    assert x > 1.0


def test_age_eligibility_proxy_uses_0_to_8_months():
    assert np.isclose(age_eligibility_proxy(1.0), 0.85)


def test_gap_formula():
    annual, gap = routine_mcv1_gap(0.5, 100_000)
    assert np.isclose(annual, 100_000 / 5 * (1 - 0.5 * 0.93))
    assert np.isclose(gap, annual * 3.5)


def test_province_seed_is_stable_and_province_specific():
    assert stable_province_seed("Uruzgan") == 2856952788
    assert stable_province_seed("Kabul") == 3031591316
    assert stable_province_seed("Uruzgan") != stable_province_seed("Kabul")
