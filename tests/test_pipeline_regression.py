from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from measles_analysis.core import DETERMINISTIC_FIGURE_PROVINCES, STOCHASTIC_PROVINCES, TABLE3_PROVINCES
from measles_analysis.pipeline import run_pipeline, validate_inputs

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"

LOCKED = {
    "direct_national_mics_mcv1_pct": 51.2,
    "direct_national_r_cov": 7.8576,
    "u5_weighted_provincial_mcv1_pct": 50.56681539402753,
    "u5_weighted_r_cov": 7.945929252533159,
    "national_u5_population": 8146035.52681,
    "routine_mcv1_unprotected_3p5yr_n": 3020631.692617184,
    "sia_preferential_n": 4056042.0399888135,
    "sia_preferential_pct_u5": 49.79160754504057,
    "sia_random_n": 7530647.856211407,
    "sia_random_pct_u5": 92.44555626386301,
    "one_dose_required_pct": 100.35842293906809,
    "two_dose_required_pct": 96.21993127147766,
    "week14_kendall_tau_b": 0.33104795088367867,
    "week14_kendall_p": 0.011157716431098409,
    "week14_spearman_rho": 0.4302580494917712,
    "week14_spearman_p": 0.011088135566639564,
    "monte_carlo_r_cov_ci_low": 7.690345305309093,
    "monte_carlo_r_cov_ci_high": 8.202344006034181,
    "monte_carlo_gap_ci_low_n": 2923471.88328927,
    "monte_carlo_gap_ci_high_n": 3118107.3315593363,
    "monte_carlo_sia_pref_ci_low_n": 3906795.0210519265,
    "monte_carlo_sia_pref_ci_high_n": 4205774.204262473,
    "monte_carlo_sia_rand_ci_low_n": 7412834.064133138,
    "monte_carlo_sia_rand_ci_high_n": 7586402.775584033,
}


def test_input_integrity():
    validate_inputs(DATA)


def test_locked_province_lists():
    assert DETERMINISTIC_FIGURE_PROVINCES == [
        "Uruzgan", "Paktika", "Helmand", "Kabul", "Herat", "Bamyan"
    ]
    assert STOCHASTIC_PROVINCES == [
        "Uruzgan", "Paktika", "Nuristan", "Herat", "Helmand",
        "Kabul", "Panjshir", "Nimroz", "Bamyan",
    ]
    assert TABLE3_PROVINCES == [
        "Uruzgan", "Paktika", "Nuristan", "Helmand", "Nimroz",
        "Panjshir", "Kabul", "Herat", "Bamyan",
    ]


def test_headline_and_monte_carlo_results(generated_outputs):
    _, key = generated_outputs
    for metric, expected in LOCKED.items():
        assert np.isclose(key[metric], expected, rtol=0, atol=1e-10), metric


def test_campaign_scenario_locked_values(generated_outputs):
    out, _ = generated_outputs
    c = pd.read_csv(out / "campaign_scenario_results.csv").set_index("SIA_reach_pct")
    expected = {
        70.0: (2.773129309134073, 1),
        90.0: (1.2951864681629044, 5),
        95.0: (0.9257007579201131, 17),
        100.0: (0.5562150476773209, 34),
    }
    for reach, (rcov, n) in expected.items():
        assert np.isclose(c.loc[reach, "U5_weighted_Rcov"], rcov, atol=1e-12)
        assert int(c.loc[reach, "provinces_Rcov_lt1"]) == n


def test_natural_history_uses_final_six_and_matches_appendix(generated_outputs):
    out, _ = generated_outputs
    n = pd.read_csv(out / "seir_natural_history_sensitivity.csv")
    assert list(dict.fromkeys(n["province"])) == DETERMINISTIC_FIGURE_PROVINCES

    expected = {
        "Uruzgan": ((26490, 31478, 24.75, 28.75), (21971, 30324, 25.00, 32.00), 92.0, 100.0),
        "Paktika": ((17093, 20318, 25.25, 29.25), (14169, 19569, 25.50, 32.50), 84.8, 100.0),
        "Helmand": ((87849, 104487, 31.50, 36.50), (72768, 100633, 31.75, 40.50), 72.6, 100.0),
        "Kabul": ((62227, 74237, 49.75, 58.00), (51305, 71456, 51.00, 64.25), 37.6, 99.6),
        "Herat": ((18190, 21746, 65.75, 77.00), (14956, 20924, 68.25, 85.25), 24.5, 97.5),
        "Bamyan": ((1012, 1212, 76.00, 89.50), (830, 1166, 79.75, 98.50), 15.9, 90.8),
    }
    for province, (lat, inf, ar_u5, ar_sus) in expected.items():
        p = n[n.province == province]
        l = p[p.vary == "latent"]
        i = p[p.vary == "infectious"]
        got_lat = (round(l.peak_infectious.min()), round(l.peak_infectious.max()), l.peak_day.min(), l.peak_day.max())
        got_inf = (round(i.peak_infectious.min()), round(i.peak_infectious.max()), i.peak_day.min(), i.peak_day.max())
        assert got_lat == lat
        assert got_inf == inf
        assert round(100 * p.attack_rate_total_u5.iloc[0], 1) == ar_u5
        assert round(100 * p.attack_rate_initial_susceptible.iloc[0], 1) == ar_sus


@pytest.mark.slow
def test_stochastic_locked_results(generated_outputs):
    out, _ = generated_outputs
    s = pd.read_csv(out / "stochastic_results_corrected.csv").set_index("province")
    expected = {
        "Uruzgan": (0.939, 26732, 26490, 26973.2),
        "Paktika": (0.918, 17257, 17061.55, 17455.75),
        "Nuristan": (0.924, 8735, 8601, 8869),
        "Herat": (0.773, 18322, 18104, 18581.6),
        "Helmand": (0.914, 88536.5, 88099.65, 88963.7),
        "Kabul": (0.832, 62633, 62196.1, 63067.35),
        "Panjshir": (0.857, 2804, 2723.6, 2886),
        "Nimroz": (0.863, 6614, 6490.1, 6745.8),
        "Bamyan": (0.659, 1046, 984, 1104),
    }
    for province, vals in expected.items():
        got = s.loc[province]
        assert np.isclose(got.outbreak_probability, vals[0], atol=0)
        assert np.isclose(got.peak_median_established, vals[1], atol=0)
        assert np.isclose(got.peak_p05_established, vals[2], atol=1e-12)
        assert np.isclose(got.peak_p95_established, vals[3], atol=1e-12)
        assert int(got.remaining_active_at_horizon) == 0
