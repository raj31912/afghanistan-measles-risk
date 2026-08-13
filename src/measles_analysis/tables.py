"""Create manuscript-ready CSV tables from exact model outputs."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .core import DETERMINISTIC_FIGURE_PROVINCES, TABLE3_PROVINCES


def _reach_category(x_pct: float) -> str:
    if x_pct >= 80:
        return "≥80%"
    if x_pct >= 60:
        return "60–<80%"
    if x_pct >= 40:
        return "40–<60%"
    return "<40%"


def _category_robust(low_pct: float, high_pct: float) -> bool:
    return _reach_category(low_pct) == _reach_category(high_pct)


def table1_display(prov: pd.DataFrame, key: dict[str, float]) -> pd.DataFrame:
    d = prov.sort_values("r_cov", ascending=False).copy()
    rows = []

    # National row: direct MICS coverage but U5-weighted model aggregate for Rcov/SIA.
    rows.append({
        "Province": "NATIONAL / aggregate",
        "Total population (2026)": f"{int(round(d.population_total_2026.sum())):,}",
        "Under-five population": f"{int(round(d.population_u5.sum())):,}",
        "MCV1 coverage % (95% CI)": "51.2 (49.1-53.2)",
        "Rcov (95% CI)": f"{key['u5_weighted_r_cov']:.2f} ({key['monte_carlo_r_cov_ci_low']:.2f}-{key['monte_carlo_r_cov_ci_high']:.2f})",
        "3.5-yr upper-bound unprotected": f"{int(round(key['routine_mcv1_unprotected_3p5yr_n']/100)*100):,}",
        "Preferential SIA reach % (95% CI)": f"{key['sia_preferential_pct_u5']:.1f} ({100*key['monte_carlo_sia_pref_ci_low_n']/key['national_u5_population']:.1f}-{100*key['monte_carlo_sia_pref_ci_high_n']/key['national_u5_population']:.1f})",
        "Random SIA reach % (95% CI)": f"{key['sia_random_pct_u5']:.1f} ({100*key['monte_carlo_sia_rand_ci_low_n']/key['national_u5_population']:.1f}-{100*key['monte_carlo_sia_rand_ci_high_n']/key['national_u5_population']:.1f})",
        "Preferential children needing SIA": f"{int(round(key['sia_preferential_n'])):,}",
        "Random children needing SIA": f"{int(round(key['sia_random_n'])):,}",
        "Preferential reach category / robust": "40–<60% / Yes",
    })

    for r in d.itertuples(index=False):
        lowp, highp = 100*r.sia_pref_ci_low, 100*r.sia_pref_ci_high
        rows.append({
            "Province": r.province,
            "Total population (2026)": f"{int(round(r.population_total_2026)):,}",
            "Under-five population": f"{int(round(r.population_u5)):,}",
            "MCV1 coverage % (95% CI)": f"{100*r.mcv1:.1f} ({100*r.mcv1_ci_low:.1f}-{100*r.mcv1_ci_high:.1f})",
            "Rcov (95% CI)": f"{r.r_cov:.2f} ({r.r_cov_ci_low:.2f}-{r.r_cov_ci_high:.2f})",
            "3.5-yr upper-bound unprotected": f"{int(round(r.routine_mcv1_unprotected_3p5yr/100)*100):,}",
            "Preferential SIA reach % (95% CI)": f"{100*r.sia_preferential:.1f} ({lowp:.1f}-{highp:.1f})",
            "Random SIA reach % (95% CI)": f"{100*r.sia_random:.1f} ({100*r.sia_rand_ci_low:.1f}-{100*r.sia_rand_ci_high:.1f})",
            "Preferential children needing SIA": f"{int(round(r.sia_preferential_n)):,}",
            "Random children needing SIA": f"{int(round(r.sia_random_n)):,}",
            "Preferential reach category / robust": f"{_reach_category(100*r.sia_preferential)} / {'Yes' if _category_robust(lowp, highp) else 'No'}",
        })
    return pd.DataFrame(rows)


def table2_display(prov: pd.DataFrame, stochastic: pd.DataFrame) -> pd.DataFrame:
    d = prov.set_index("province").loc[DETERMINISTIC_FIGURE_PROVINCES].reset_index()
    s = stochastic.set_index("province")
    rows=[]
    for r in d.itertuples(index=False):
        sr=s.loc[r.province]
        rows.append({
            "Province": r.province,
            "Under-five pop. (N)": f"{int(round(r.population_u5)):,}",
            "MCV1 (%)": f"{100*r.mcv1:.1f}%",
            "Rcov": f"{r.r_cov:.2f}",
            "Deterministic peak I": f"{r.peak_infectious:,.0f}",
            "Peak day": f"{r.peak_day:.1f}",
            "Attack rate (total U5)": f"{100*r.attack_rate_total_u5:.1f}%",
            "Attack rate (initial susceptible)": f"{100*r.attack_rate_initial_susceptible:.1f}%",
            "Stochastic peak median": f"{sr.peak_median_established:.0f}",
            "P05-P95 (stochastic)": f"{sr.peak_p05_established:.0f}-{sr.peak_p95_established:.0f}",
            "Establishment probability (%)": f"{100*sr.outbreak_probability:.1f}%",
        })
    return pd.DataFrame(rows)


def table3_display(prov: pd.DataFrame, stochastic: pd.DataFrame) -> pd.DataFrame:
    d=prov.set_index("province")
    rows=[]
    for sr in stochastic.set_index("province").loc[TABLE3_PROVINCES].itertuples():
        r=d.loc[sr.Index]
        rows.append({
            "Province": sr.Index,
            "MCV1 (%)": f"{100*r.mcv1:.1f}%",
            "Rcov": f"{r.r_cov:.2f}",
            "Under-five pop.": f"{int(round(r.population_u5)):,}",
            "Deterministic peak I": f"{r.peak_infectious:,.0f}",
            "Stochastic peak median": f"{sr.peak_median_established:.0f}",
            "P05-P95": f"{sr.peak_p05_established:.0f}-{sr.peak_p95_established:.0f}",
            "Establishment probability (%)": f"{100*sr.outbreak_probability:.1f}",
            "Non-establishment probability (%)": f"{100*sr.non_establishment_probability:.1f}",
            "Median attack rate (total U5)": f"{100*sr.attack_total_u5_median_established:.1f}%",
        })
    return pd.DataFrame(rows)
