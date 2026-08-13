"""Core model functions for the Afghanistan measles Vaccine revision.

All calculations in this module are deterministic given their inputs and documented
random seeds. The functions intentionally distinguish:

* direct national survey estimates from U5-weighted provincial aggregates;
* survey sampling uncertainty from demographic stochasticity;
* an idealized preferential-targeting SIA lower bound from random U5 targeting;
* the coverage-adjusted transmission-potential index Rcov from an empirically
  estimated province-specific effective reproduction number.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import truncnorm


@dataclass(frozen=True)
class Parameters:
    r0: float = 15.0
    r0_low: float = 12.0
    r0_high: float = 18.0
    ve_one_dose: float = 0.93
    ve_two_dose: float = 0.97
    latent_days: float = 8.0
    infectious_days: float = 4.0
    horizon_years_gap: float = 3.5
    seir_days: float = 365.0
    seir_step_days: float = 0.25
    stochastic_sims: int = 1000
    stochastic_dt: float = 0.1
    stochastic_seed: int = 42
    outbreak_threshold: int = 10
    monte_carlo_draws: int = 100_000
    monte_carlo_seed: int = 20260810

    @property
    def hit(self) -> float:
        return 1.0 - 1.0 / self.r0

    @property
    def sigma(self) -> float:
        return 1.0 / self.latent_days

    @property
    def gamma(self) -> float:
        return 1.0 / self.infectious_days

    @property
    def beta(self) -> float:
        return self.r0 * self.gamma


P = Parameters()


def coverage_adjusted_transmission(coverage, r0: float = P.r0, ve: float = P.ve_one_dose):
    """Coverage-adjusted transmission-potential index, Rcov.

    This is a model-derived coverage index under homogeneous mixing, not an
    empirically estimated province-specific effective reproduction number.
    """
    coverage = np.asarray(coverage, dtype=float)
    return r0 * (1.0 - coverage * ve)


def preferential_sia(coverage, r0: float = P.r0, ve: float = P.ve_one_dose):
    """Idealized preferential-targeting SIA reach requirement.

    Doses are assumed to reach children lacking MCV1-derived protection first.
    Values are deliberately not capped at 100%; >100% indicates one-pass
    infeasibility under the assumptions.
    """
    coverage = np.asarray(coverage, dtype=float)
    p_immune = coverage * ve
    hit = 1.0 - 1.0 / r0
    return np.maximum((hit - p_immune) / ve, 0.0)


def random_sia(coverage, r0: float = P.r0, ve: float = P.ve_one_dose):
    """Random U5-targeting SIA reach requirement.

    Reach is random within the modeled under-five population irrespective of
    existing MCV1-derived protection. Values are not capped at 100%.
    """
    coverage = np.asarray(coverage, dtype=float)
    p_immune = coverage * ve
    hit = 1.0 - 1.0 / r0
    denom = (1.0 - p_immune) * ve
    raw = np.divide(
        hit - p_immune,
        denom,
        out=np.zeros_like(p_immune, dtype=float),
        where=denom > 0,
    )
    return np.maximum(raw, 0.0)


def random_sia_post_campaign_rcov(coverage, campaign_reach, r0: float = P.r0, ve: float = P.ve_one_dose):
    """Post-campaign Rcov under random SIA reach within U5."""
    coverage = np.asarray(coverage, dtype=float)
    p_current = coverage * ve
    p_post = p_current + (1.0 - p_current) * float(campaign_reach) * ve
    return r0 * (1.0 - p_post)


def age_eligibility_proxy(coverage, ineligible_months: int = 9):
    """Simple U5 proxy assigning zero MCV1-derived immunity below eligibility.

    Routine MCV1 begins at approximately 9 months. With a uniform 0--59 month
    U5 age distribution, children aged 0--8 months are 15% of U5. Maternal
    antibody is not represented, so this is a sensitivity proxy rather than an
    age-structured transmission model.
    """
    if not 0 <= ineligible_months <= 60:
        raise ValueError("ineligible_months must be between 0 and 60")
    return np.asarray(coverage, dtype=float) * (1.0 - ineligible_months / 60.0)


def routine_mcv1_gap(coverage, population_u5, years: float = P.horizon_years_gap, ve: float = P.ve_one_dose):
    """Standardized upper-bound routine-MCV1 protection-gap approximation."""
    coverage = np.asarray(coverage, dtype=float)
    pop = np.asarray(population_u5, dtype=float)
    annual_cohort = pop / 5.0
    annual_unprotected = annual_cohort * (1.0 - coverage * ve)
    return annual_unprotected, annual_unprotected * years


def solve_seir(coverage: float, population: float, *, params: Parameters = P) -> dict[str, Any]:
    """Solve the deterministic U5 SEIR ODE system forward in time.

    One imported infectious case is removed from the susceptible pool so mass
    is conserved. Incident infections caused after introduction are measured by
    S(0)-S(T), excluding the imported seed case.
    """
    N = float(population)
    immune_fraction = float(np.clip(coverage * params.ve_one_dose, 0.0, 1.0))
    susceptible_before_intro = N * (1.0 - immune_fraction)
    I0 = 1.0
    if susceptible_before_intro <= I0:
        raise ValueError("Population has too few susceptible individuals for one imported case.")

    S0 = susceptible_before_intro - I0
    E0 = 0.0
    R0_init = N - susceptible_before_intro

    def rhs(_t, y):
        S, E, I, R = y
        new_exp = params.beta * S * I / N
        return [
            -new_exp,
            new_exp - params.sigma * E,
            params.sigma * E - params.gamma * I,
            params.gamma * I,
        ]

    t_eval = np.arange(0.0, params.seir_days + params.seir_step_days / 2, params.seir_step_days)
    sol = solve_ivp(
        rhs,
        (0.0, params.seir_days),
        [S0, E0, I0, R0_init],
        t_eval=t_eval,
        rtol=1e-8,
        atol=1e-8,
        method="RK45",
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    S, E, I, R = sol.y
    if not np.allclose(S + E + I + R, N, rtol=1e-8, atol=1e-3):
        raise RuntimeError("Deterministic SEIR mass-balance check failed.")

    peak_idx = int(np.argmax(I))
    incident = float(S0 - S[-1])
    attack_u5 = incident / N
    attack_sus = incident / S0
    if not (0 <= attack_u5 <= 1.000001 and 0 <= attack_sus <= 1.000001):
        raise RuntimeError("Deterministic attack-rate bounds check failed.")

    return {
        "time": sol.t,
        "S": S,
        "E": E,
        "I": I,
        "R": R,
        "summary": {
            "peak_infectious": float(I[peak_idx]),
            "peak_day": float(sol.t[peak_idx]),
            "incident_infections_after_introduction": incident,
            "attack_rate_total_u5": attack_u5,
            "attack_rate_initial_susceptible": attack_sus,
        },
    }


def run_seir_summary(coverage: float, population: float, *, params: Parameters = P) -> dict[str, float]:
    return solve_seir(coverage, population, params=params)["summary"]


STOCHASTIC_PROVINCES = [
    "Uruzgan", "Paktika", "Nuristan", "Herat", "Helmand",
    "Kabul", "Panjshir", "Nimroz", "Bamyan",
]

DETERMINISTIC_FIGURE_PROVINCES = [
    "Uruzgan", "Paktika", "Helmand", "Kabul", "Herat", "Bamyan",
]

# Ordering used in the final manuscript Table 3. This is deliberately separate
# from STOCHASTIC_PROVINCES, whose order is retained for locked simulation output.
TABLE3_PROVINCES = [
    "Uruzgan", "Paktika", "Nuristan", "Helmand", "Nimroz",
    "Panjshir", "Kabul", "Herat", "Bamyan",
]


def stable_province_seed(province: str, base_seed: int = P.stochastic_seed) -> int:
    """Create deterministic, province-specific RNG streams from one base seed."""
    digest = hashlib.sha256(f"{base_seed}:{province}".encode()).digest()
    return int.from_bytes(digest[:8], "little") % (2**32)


def run_stochastic_summary(coverage: float, population: float, province: str, *, params: Parameters = P) -> dict[str, float]:
    """Bounded-Poisson tau-leaping summary for one province.

    Establishment is operationally defined as peak I > outbreak_threshold.
    Peak percentiles are conditional on establishment and quantify demographic
    stochasticity, not survey sampling uncertainty.
    """
    rng = np.random.default_rng(stable_province_seed(province, params.stochastic_seed))
    N = int(round(population))
    immune_fraction = float(np.clip(coverage * params.ve_one_dose, 0.0, 1.0))
    susceptible_before_intro = int(round(N * (1.0 - immune_fraction)))
    if susceptible_before_intro <= 1:
        raise ValueError("Population has too few susceptible individuals for one imported case.")

    S0 = susceptible_before_intro - 1
    S = np.full(params.stochastic_sims, S0, dtype=np.int64)
    E = np.zeros(params.stochastic_sims, dtype=np.int64)
    I = np.ones(params.stochastic_sims, dtype=np.int64)
    R = np.full(params.stochastic_sims, N - susceptible_before_intro, dtype=np.int64)
    peak = I.copy()
    peak_step = np.zeros(params.stochastic_sims, dtype=np.int32)

    n_steps = int(round(params.seir_days / params.stochastic_dt))
    for step in range(1, n_steps + 1):
        active = (I > 0) | (E > 0)
        if not np.any(active):
            break
        idx = np.flatnonzero(active)
        Sa, Ea, Ia = S[idx], E[idx], I[idx]
        exposures = np.minimum(rng.poisson(params.beta * Sa * Ia / N * params.stochastic_dt), Sa)
        progressions = np.minimum(rng.poisson(params.sigma * Ea * params.stochastic_dt), Ea)
        recoveries = np.minimum(rng.poisson(params.gamma * Ia * params.stochastic_dt), Ia)

        S[idx] = Sa - exposures
        E[idx] = Ea + exposures - progressions
        I[idx] = Ia + progressions - recoveries
        R[idx] = R[idx] + recoveries

        higher = I[idx] > peak[idx]
        if np.any(higher):
            hidx = idx[higher]
            peak[hidx] = I[hidx]
            peak_step[hidx] = step

    if not np.all(S + E + I + R == N):
        raise RuntimeError(f"Stochastic SEIR mass-balance check failed for {province}.")

    established = peak > params.outbreak_threshold
    established_peaks = peak[established]
    established_peak_days = peak_step[established] * params.stochastic_dt
    incident = S0 - S
    attack_u5 = incident / N
    attack_sus = incident / S0

    def q(arr, prob):
        return float(np.quantile(arr, prob)) if len(arr) else np.nan

    return {
        "province": province,
        "province_seed": stable_province_seed(province, params.stochastic_seed),
        "n_sims": params.stochastic_sims,
        "outbreak_probability": float(established.mean()),
        "non_establishment_probability": float(1.0 - established.mean()),
        "n_established": int(established.sum()),
        "peak_median_established": q(established_peaks, 0.50),
        "peak_p05_established": q(established_peaks, 0.05),
        "peak_p95_established": q(established_peaks, 0.95),
        "peak_day_median_established": q(established_peak_days, 0.50),
        "attack_total_u5_median_established": q(attack_u5[established], 0.50),
        "attack_initial_sus_median_established": q(attack_sus[established], 0.50),
        "remaining_active_at_horizon": int(np.sum((I > 0) | (E > 0))),
    }


def draw_truncated_normal(mean: float, se: float, n: int, rng) -> np.ndarray:
    if se <= 0:
        return np.full(n, mean, dtype=float)
    a = (0.0 - mean) / se
    b = (1.0 - mean) / se
    return truncnorm.rvs(a, b, loc=mean, scale=se, size=n, random_state=rng)
