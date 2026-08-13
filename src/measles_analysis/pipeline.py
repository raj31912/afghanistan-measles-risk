"""End-to-end audited analysis pipeline for the Afghanistan measles revision."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

from .core import (
    P,
    Parameters,
    STOCHASTIC_PROVINCES,
    DETERMINISTIC_FIGURE_PROVINCES,
    age_eligibility_proxy,
    coverage_adjusted_transmission,
    draw_truncated_normal,
    preferential_sia,
    random_sia,
    random_sia_post_campaign_rcov,
    routine_mcv1_gap,
    run_seir_summary,
    run_stochastic_summary,
    solve_seir,
)
from .tables import table1_display, table2_display, table3_display


REQUIRED_INPUTS = {
    "mics_official_sampling_errors.csv": {
        "province", "mcv1", "mcv1_se", "mcv1_ci_low", "mcv1_ci_high"
    },
    "mics_national_official.csv": {
        "mcv1", "mcv1_se", "mcv1_ci_low", "mcv1_ci_high", "weighted_n", "unweighted_n"
    },
    "population_2026_provincial_clean.csv": {
        "province", "population_u5", "population_total_2026"
    },
    "who_week14_measles_incidence_categories.csv": {
        "province", "incidence_category"
    },
}

EXPECTED_PROVINCES = {
    "Badakhshan", "Badghis", "Baghlan", "Balkh", "Bamyan", "Daykundi",
    "Farah", "Faryab", "Ghazni", "Ghor", "Helmand", "Herat", "Jowzjan",
    "Kabul", "Kandahar", "Kapisa", "Khost", "Kunar", "Kunduz", "Laghman",
    "Logar", "Nangarhar", "Nimroz", "Nuristan", "Paktia", "Paktika",
    "Panjshir", "Parwan", "Samangan", "Sar-e Pul", "Takhar", "Uruzgan",
    "Wardak", "Zabul",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_inputs(data_dir: Path) -> None:
    for fname, required_cols in REQUIRED_INPUTS.items():
        path = data_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"Required input missing: {path}")
        df = pd.read_csv(path)
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"{fname} is missing required columns: {sorted(missing)}")

    mics = pd.read_csv(data_dir / "mics_official_sampling_errors.csv")
    pop = pd.read_csv(data_dir / "population_2026_provincial_clean.csv")
    who = pd.read_csv(data_dir / "who_week14_measles_incidence_categories.csv")

    for label, df in [("MICS", mics), ("population", pop), ("WHO", who)]:
        if len(df) != 34 or df["province"].nunique() != 34:
            raise ValueError(f"{label} input must contain 34 unique provinces; found {len(df)} rows / {df['province'].nunique()} unique.")
        names = set(df["province"])
        if names != EXPECTED_PROVINCES:
            raise ValueError(f"{label} province set differs from the locked 34-province set. Missing={sorted(EXPECTED_PROVINCES-names)}; extra={sorted(names-EXPECTED_PROVINCES)}")

    if not ((mics["mcv1"] >= 0) & (mics["mcv1"] <= 1)).all():
        raise ValueError("MCV1 coverage must be within [0,1].")
    if not ((mics["mcv1_ci_low"] >= 0) & (mics["mcv1_ci_high"] <= 1) & (mics["mcv1_ci_low"] <= mics["mcv1"]) & (mics["mcv1"] <= mics["mcv1_ci_high"])).all():
        raise ValueError("MICS confidence intervals are invalid.")
    if not (mics["mcv1_se"] > 0).all():
        raise ValueError("MICS provincial SE values must be positive.")
    if not (pop["population_u5"] > 0).all():
        raise ValueError("U5 populations must be positive.")
    if not (who["incidence_category"].isin([1,2,3,4,5])).all():
        raise ValueError("WHO incidence categories must be ordinal integers 1-5.")


def load_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.Series]:
    validate_inputs(data_dir)
    mics = pd.read_csv(data_dir / "mics_official_sampling_errors.csv")
    national = pd.read_csv(data_dir / "mics_national_official.csv").iloc[0]
    pop = pd.read_csv(data_dir / "population_2026_provincial_clean.csv")
    who = pd.read_csv(data_dir / "who_week14_measles_incidence_categories.csv")

    # Preserve the MICS input order. This matters only for exact reproducibility of
    # seeded Monte Carlo draws; analytical estimates are invariant to row order.
    data = mics.merge(pop, on="province", validate="one_to_one", sort=False).merge(
        who, on="province", validate="one_to_one", sort=False
    )
    if len(data) != 34:
        raise RuntimeError("Merged analysis dataset did not retain 34 provinces.")
    return data, national


def compute_provincial_results(data: pd.DataFrame, params: Parameters = P) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = data.copy()
    d["immune_fraction_mcv1"] = d["mcv1"] * params.ve_one_dose
    d["r_cov"] = coverage_adjusted_transmission(d["mcv1"], params.r0, params.ve_one_dose)
    d["r_cov_r0_12"] = coverage_adjusted_transmission(d["mcv1"], params.r0_low, params.ve_one_dose)
    d["r_cov_r0_18"] = coverage_adjusted_transmission(d["mcv1"], params.r0_high, params.ve_one_dose)
    d["sia_preferential"] = preferential_sia(d["mcv1"], params.r0, params.ve_one_dose)
    d["sia_random"] = random_sia(d["mcv1"], params.r0, params.ve_one_dose)
    d["sia_preferential_n"] = d["sia_preferential"] * d["population_u5"]
    d["sia_random_n"] = d["sia_random"] * d["population_u5"]
    d["sia_preferential_unattainable_one_pass"] = d["sia_preferential"] > 1.0
    d["sia_random_unattainable_one_pass"] = d["sia_random"] > 1.0

    d["mcv1_age_eligibility_proxy"] = age_eligibility_proxy(d["mcv1"], ineligible_months=9)
    d["r_cov_age_eligibility_proxy"] = coverage_adjusted_transmission(
        d["mcv1_age_eligibility_proxy"], params.r0, params.ve_one_dose
    )

    annual_gap, gap_35 = routine_mcv1_gap(
        d["mcv1"], d["population_u5"], params.horizon_years_gap, params.ve_one_dose
    )
    d["annual_routine_mcv1_unprotected"] = annual_gap
    d["routine_mcv1_unprotected_3p5yr"] = gap_35

    # Province-level design-based uncertainty propagated through monotone functions.
    d["r_cov_ci_low"] = coverage_adjusted_transmission(d["mcv1_ci_high"], params.r0, params.ve_one_dose)
    d["r_cov_ci_high"] = coverage_adjusted_transmission(d["mcv1_ci_low"], params.r0, params.ve_one_dose)
    d["sia_pref_ci_low"] = preferential_sia(d["mcv1_ci_high"], params.r0, params.ve_one_dose)
    d["sia_pref_ci_high"] = preferential_sia(d["mcv1_ci_low"], params.r0, params.ve_one_dose)
    d["sia_rand_ci_low"] = random_sia(d["mcv1_ci_high"], params.r0, params.ve_one_dose)
    d["sia_rand_ci_high"] = random_sia(d["mcv1_ci_low"], params.r0, params.ve_one_dose)
    # The protection gap decreases monotonically with MCV1 coverage, so the
    # lower gap bound uses the upper coverage CI and vice versa.
    _, d["routine_gap_ci_low"] = routine_mcv1_gap(
        d["mcv1_ci_high"], d["population_u5"], params.horizon_years_gap, params.ve_one_dose
    )
    _, d["routine_gap_ci_high"] = routine_mcv1_gap(
        d["mcv1_ci_low"], d["population_u5"], params.horizon_years_gap, params.ve_one_dose
    )

    # Deterministic SEIR summary for every province.
    seir_rows=[]
    for row in d.itertuples(index=False):
        seir_rows.append({"province": row.province, **run_seir_summary(row.mcv1, row.population_u5, params=params)})
    seir = pd.DataFrame(seir_rows)
    d = d.merge(seir, on="province", validate="one_to_one", sort=False)
    return d, seir


def compute_monte_carlo(data: pd.DataFrame, params: Parameters = P) -> dict[str, float]:
    rng = np.random.default_rng(params.monte_carlo_seed)
    draws = np.column_stack([
        draw_truncated_normal(r.mcv1, r.mcv1_se, params.monte_carlo_draws, rng)
        for r in data.itertuples(index=False)
    ])
    weights = data["population_u5"].to_numpy(float)
    weighted_cov_draw = np.average(draws, axis=1, weights=weights)
    rcov_draw = coverage_adjusted_transmission(weighted_cov_draw, params.r0, params.ve_one_dose)
    pref_draw = (preferential_sia(draws, params.r0, params.ve_one_dose) * weights).sum(axis=1)
    rand_draw = (random_sia(draws, params.r0, params.ve_one_dose) * weights).sum(axis=1)
    _, gap_matrix = routine_mcv1_gap(draws, weights, params.horizon_years_gap, params.ve_one_dose)
    gap_draw = gap_matrix.sum(axis=1)

    def ci(x):
        lo, hi = np.quantile(x, [0.025,0.975])
        return float(lo), float(hi)

    cov_lo,cov_hi=ci(weighted_cov_draw)
    rc_lo,rc_hi=ci(rcov_draw)
    gap_lo,gap_hi=ci(gap_draw)
    pref_lo,pref_hi=ci(pref_draw)
    rand_lo,rand_hi=ci(rand_draw)
    return {
        "monte_carlo_u5_weighted_mcv1_ci_low_pct": cov_lo*100,
        "monte_carlo_u5_weighted_mcv1_ci_high_pct": cov_hi*100,
        "monte_carlo_r_cov_ci_low": rc_lo,
        "monte_carlo_r_cov_ci_high": rc_hi,
        "monte_carlo_gap_ci_low_n": gap_lo,
        "monte_carlo_gap_ci_high_n": gap_hi,
        "monte_carlo_sia_pref_ci_low_n": pref_lo,
        "monte_carlo_sia_pref_ci_high_n": pref_hi,
        "monte_carlo_sia_rand_ci_low_n": rand_lo,
        "monte_carlo_sia_rand_ci_high_n": rand_hi,
    }


def compute_structural_sensitivity(data: pd.DataFrame, params: Parameters = P) -> pd.DataFrame:
    weights=data["population_u5"].to_numpy(float)
    total_u5=float(weights.sum())
    rows=[]
    specs=[
        ("Reference", params.r0, params.ve_one_dose, 1.0),
        ("Lower R0", params.r0_low, params.ve_one_dose, 1.0),
        ("Upper R0", params.r0_high, params.ve_one_dose, 1.0),
        ("Higher VE1", params.r0, 0.95, 1.0),
        ("Age-eligibility proxy (0–8 months zero MCV1-derived immunity)", params.r0, params.ve_one_dose, 0.85),
        ("Best-favourable R0+VE1", params.r0_low, 0.95, 1.0),
    ]
    for label,r0_s,ve1_s,age_factor in specs:
        cov=data["mcv1"].to_numpy(float)*age_factor
        rc=coverage_adjusted_transmission(cov,r0_s,ve1_s)
        pref=preferential_sia(cov,r0_s,ve1_s)
        rand=random_sia(cov,r0_s,ve1_s)
        rows.append({
            "scenario":label,"r0":r0_s,"ve1":ve1_s,"age_factor":age_factor,
            "u5_weighted_mcv1_pct":np.average(cov,weights=weights)*100,
            "u5_weighted_r_cov":np.average(rc,weights=weights),
            "provinces_r_cov_gt_1":int(np.sum(rc>1)),
            "one_dose_required_pct":(1-1/r0_s)/ve1_s*100,
            "two_dose_required_pct_ve2_0p97":(1-1/r0_s)/params.ve_two_dose*100,
            "sia_preferential_n":float(np.sum(pref*weights)),
            "sia_preferential_pct_u5":float(np.sum(pref*weights)/total_u5*100),
            "sia_random_n":float(np.sum(rand*weights)),
            "sia_random_pct_u5":float(np.sum(rand*weights)/total_u5*100),
            "n_provinces_random_requirement_gt_100pct":int(np.sum(rand>1)),
        })
    return pd.DataFrame(rows)


def compute_two_dose_sensitivity(params: Parameters = P) -> pd.DataFrame:
    rows=[]
    for r0_s in (params.r0_low,params.r0,params.r0_high):
        for ve2_s in (0.96,0.97,0.99):
            rows.append({"r0":r0_s,"ve2":ve2_s,"two_dose_required_pct":(1-1/r0_s)/ve2_s*100})
    return pd.DataFrame(rows)


def compute_natural_history_sensitivity(data: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity for the same six provinces used in final Table 2/Figure 3."""
    by=data.set_index("province")
    rows=[]
    for prov in DETERMINISTIC_FIGURE_PROVINCES:
        cov=float(by.loc[prov,"mcv1"]); pop=float(by.loc[prov,"population_u5"])
        for latent in (6.0,7.0,8.0):
            p=Parameters(latent_days=latent,infectious_days=4.0)
            rows.append({"province":prov,"vary":"latent","latent_days":latent,"infectious_days":4.0,**run_seir_summary(cov,pop,params=p)})
        for infectious in (3.0,4.0,5.0):
            p=Parameters(latent_days=8.0,infectious_days=infectious)
            rows.append({"province":prov,"vary":"infectious","latent_days":8.0,"infectious_days":infectious,**run_seir_summary(cov,pop,params=p)})
    return pd.DataFrame(rows)


def compute_campaign_scenario(data: pd.DataFrame, params: Parameters = P) -> pd.DataFrame:
    reaches=np.array([0,20,40,50,60,70,80,90,95,100],dtype=float)
    weights=data["population_u5"].to_numpy(float)
    coverage=data["mcv1"].to_numpy(float)
    rows=[]
    for pct in reaches:
        rc=random_sia_post_campaign_rcov(coverage,pct/100.0,params.r0,params.ve_one_dose)
        rows.append({
            "SIA_reach_pct":pct,
            "U5_weighted_Rcov":float(np.average(rc,weights=weights)),
            "provinces_Rcov_lt1":int(np.sum(rc<1)),
        })
    return pd.DataFrame(rows)


def compute_stochastic(data: pd.DataFrame, params: Parameters = P) -> pd.DataFrame:
    by=data.set_index("province")
    rows=[]
    for prov in STOCHASTIC_PROVINCES:
        row=by.loc[prov]
        s=run_stochastic_summary(float(row.mcv1),float(row.population_u5),prov,params=params)
        s.update({
            "mcv1":float(row.mcv1),"population_u5":float(row.population_u5),"r_cov":float(row.r_cov),
            "deterministic_peak":float(row.peak_infectious),"deterministic_peak_day":float(row.peak_day),
            "deterministic_attack_total_u5":float(row.attack_rate_total_u5),
        })
        rows.append(s)
    return pd.DataFrame(rows)


def write_deterministic_trajectories(data: pd.DataFrame, output_dir: Path, params: Parameters = P) -> None:
    by=data.set_index("province")
    rows=[]
    for prov in DETERMINISTIC_FIGURE_PROVINCES:
        row=by.loc[prov]
        sol=solve_seir(float(row.mcv1),float(row.population_u5),params=params)
        for t,S,E,I,R in zip(sol["time"],sol["S"],sol["E"],sol["I"],sol["R"]):
            rows.append({"province":prov,"day":t,"S":S,"E":E,"I":I,"R":R})
    pd.DataFrame(rows).to_csv(output_dir/"figure3_seir_trajectories.csv",index=False)


def run_pipeline(data_dir: Path, output_dir: Path, *, run_stochastic: bool = True, params: Parameters = P) -> dict[str, float]:
    output_dir.mkdir(parents=True,exist_ok=True)
    data,national=load_inputs(data_dir)
    prov,seir=compute_provincial_results(data,params=params)

    tau=kendalltau(prov["r_cov"],prov["incidence_category"],variant="b")
    rho=spearmanr(prov["r_cov"],prov["incidence_category"])

    direct_cov=float(national.mcv1)
    direct_rcov=float(coverage_adjusted_transmission(direct_cov,params.r0,params.ve_one_dose))
    total_u5=float(prov["population_u5"].sum())
    weighted_cov=float(np.average(prov["mcv1"],weights=prov["population_u5"]))
    weighted_rcov=float(coverage_adjusted_transmission(weighted_cov,params.r0,params.ve_one_dose))
    pref_n=float(prov["sia_preferential_n"].sum())
    rand_n=float(prov["sia_random_n"].sum())
    gap_n=float(prov["routine_mcv1_unprotected_3p5yr"].sum())

    key={
        "direct_national_mics_mcv1_pct":direct_cov*100,
        "direct_national_r_cov":direct_rcov,
        "u5_weighted_provincial_mcv1_pct":weighted_cov*100,
        "u5_weighted_r_cov":weighted_rcov,
        "national_u5_population":total_u5,
        "routine_mcv1_unprotected_3p5yr_n":gap_n,
        "sia_preferential_n":pref_n,
        "sia_preferential_pct_u5":pref_n/total_u5*100,
        "sia_random_n":rand_n,
        "sia_random_pct_u5":rand_n/total_u5*100,
        "one_dose_required_pct":params.hit/params.ve_one_dose*100,
        "two_dose_required_pct":params.hit/params.ve_two_dose*100,
        "week14_kendall_tau_b":float(tau.statistic),
        "week14_kendall_p":float(tau.pvalue),
        "week14_spearman_rho":float(rho.statistic),
        "week14_spearman_p":float(rho.pvalue),
    }
    key.update(compute_monte_carlo(prov,params=params))

    structural=compute_structural_sensitivity(prov,params=params)
    ve2=compute_two_dose_sensitivity(params=params)
    natural=compute_natural_history_sensitivity(prov)
    campaign=compute_campaign_scenario(prov,params=params)

    if run_stochastic:
        stochastic=compute_stochastic(prov,params=params)
    else:
        stochastic=pd.DataFrame()

    # Exact outputs.
    prov.to_csv(output_dir/"provincial_results_recomputed.csv",index=False)
    seir.to_csv(output_dir/"deterministic_seir_results.csv",index=False)
    structural.to_csv(output_dir/"structural_sensitivity_results.csv",index=False)
    ve2.to_csv(output_dir/"two_dose_ve_sensitivity.csv",index=False)
    natural.to_csv(output_dir/"seir_natural_history_sensitivity.csv",index=False)
    campaign.to_csv(output_dir/"campaign_scenario_results.csv",index=False)
    pd.DataFrame(key.items(),columns=["metric","value"]).to_csv(output_dir/"key_results_recomputed.csv",index=False)
    write_deterministic_trajectories(prov,output_dir,params=params)

    if run_stochastic:
        stochastic.to_csv(output_dir/"stochastic_results_corrected.csv",index=False)
        table1_display(prov,key).to_csv(output_dir/"table1_manuscript.csv",index=False)
        table2_display(prov,stochastic).to_csv(output_dir/"table2_manuscript.csv",index=False)
        table3_display(prov,stochastic).to_csv(output_dir/"table3_manuscript.csv",index=False)

    # Machine-readable reproducibility metadata.
    input_hashes={name:sha256_file(data_dir/name) for name in REQUIRED_INPUTS}
    metadata={
        "analysis_parameters":asdict(params),
        "input_sha256":input_hashes,
        "province_count":34,
        "stochastic_run":bool(run_stochastic),
        "important_interpretation":{
            "Rcov":"coverage-adjusted transmission-potential index; not empirical province-specific Re",
            "routine_gap":"standardized upper-bound routine-MCV1 protection gap; not actual susceptible population",
            "preferential_SIA":"idealized lower-bound targeting assumption",
            "stochastic_intervals":"demographic stochasticity conditional on fixed inputs; not MICS confidence intervals",
            "WHO_comparison":"retrospective ordinal geographic concordance; not external validation",
        },
    }
    with (output_dir/"analysis_metadata.json").open("w",encoding="utf-8") as f:
        json.dump(metadata,f,indent=2,sort_keys=True)

    return key
