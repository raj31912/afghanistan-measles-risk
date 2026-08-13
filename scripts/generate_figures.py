#!/usr/bin/env python3
"""Regenerate manuscript figures from audited analysis outputs.

Non-map figures (1–7) require only the repository's generated CSVs. Map figures
(7A–7C) additionally require a GADM v4.1 Afghanistan level-1 boundary file,
which is not redistributed by this repository. Pass it with --gadm.

The purpose of this script is scientific/analytical reproducibility. Journal
layout may subsequently apply cosmetic resizing without changing plotted data.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import unicodedata

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from measles_analysis.core import (
    P, DETERMINISTIC_FIGURE_PROVINCES, TABLE3_PROVINCES
)


def _save(fig, out: Path, name: str):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def _load(tables: Path):
    p = pd.read_csv(tables / "provincial_results_recomputed.csv")
    k = pd.read_csv(tables / "key_results_recomputed.csv").set_index("metric")["value"]
    s = pd.read_csv(tables / "stochastic_results_corrected.csv")
    c = pd.read_csv(tables / "campaign_scenario_results.csv")
    tr = pd.read_csv(tables / "figure3_seir_trajectories.csv")
    return p, k, s, c, tr


def figure1(p, k, out):
    d = p.sort_values("r_cov", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 11))
    xerr = np.vstack([d.r_cov - d.r_cov_ci_low, d.r_cov_ci_high - d.r_cov])
    ax.errorbar(d.r_cov, np.arange(len(d)), xerr=xerr, fmt="o", capsize=2)
    ax.set_yticks(np.arange(len(d)), d.province)
    ax.axvline(1, linestyle="--", linewidth=1)
    ax.set_xlabel("Coverage-adjusted transmission potential (Rcov)")
    ax.set_title("Coverage-adjusted Transmission Potential by Province\nAfghanistan, 2022–23 MICS")
    ax.grid(axis="x", alpha=.25)
    _save(fig, out, "Figure1_Rcov_by_province")


def figure2(p, k, out):
    d = p.sort_values("mcv1", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 11))
    y = np.arange(len(d))
    err = np.vstack([100*(d.mcv1-d.mcv1_ci_low), 100*(d.mcv1_ci_high-d.mcv1)])
    ax.barh(y, 100*d.mcv1)
    ax.errorbar(100*d.mcv1, y, xerr=err, fmt="none", capsize=2)
    ax.set_yticks(y, d.province)
    ax.axvline(float(k["direct_national_mics_mcv1_pct"]), linestyle="--", linewidth=1, label="Direct national MICS 51.2%")
    ax.set_xlim(0, 105)
    ax.set_xlabel("MCV1 coverage among children aged 12–23 months (%)")
    ax.set_title("Provincial MCV1 Coverage — Afghanistan 2022–23 MICS")
    ax.legend()
    ax.grid(axis="x", alpha=.25)
    _save(fig, out, "Figure2_MCV1_coverage")


def figure3(tr, out):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.8))
    for prov in DETERMINISTIC_FIGURE_PROVINCES:
        d = tr[tr.province == prov]
        ax1.plot(d.day, d.I, label=prov)
        ax2.plot(d.day, d.S, label=prov)
    ax1.set_title("Infectious compartment")
    ax2.set_title("Susceptible compartment")
    for ax in (ax1, ax2):
        ax.set_xlabel("Days after one imported infectious case")
        ax.grid(alpha=.2)
    ax1.set_ylabel("Children")
    ax2.set_ylabel("Children")
    ax1.legend(fontsize=8)
    ax2.legend(fontsize=8)
    fig.suptitle("Deterministic SEIR Trajectories — Six Representative Provinces")
    _save(fig, out, "Figure3_SEIR_trajectories")


def figure4(p, out):
    d = p.sort_values("routine_mcv1_unprotected_3p5yr", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 11))
    y=np.arange(len(d))
    x=d.routine_mcv1_unprotected_3p5yr
    err=np.vstack([x-d.routine_gap_ci_low, d.routine_gap_ci_high-x])
    ax.barh(y, x/1000)
    ax.errorbar(x/1000, y, xerr=err/1000, fmt="none", capsize=2)
    ax.set_yticks(y, d.province)
    ax.set_xlabel("3.5-year upper-bound routine-MCV1 unprotected accumulation (thousands)")
    ax.set_title("Provincial Contribution to the Standardized Routine-MCV1 Protection Gap")
    ax.grid(axis="x", alpha=.25)
    _save(fig, out, "Figure4_upper_bound_routine_gap")


def figure5(p, out):
    d = p.sort_values("sia_random", ascending=True)
    y=np.arange(len(d))
    fig, ax = plt.subplots(figsize=(11, 12))
    ax.barh(y, 100*d.sia_random, label="Random targeting")
    ax.scatter(100*d.sia_preferential, y, marker="|", s=100, label="Preferential lower bound")
    ax.set_yticks(y, d.province)
    ax.set_xlabel("Required SIA reach (% of modeled U5 population)")
    ax.set_title("Minimum SIA Reach Required to Reduce Rcov Below 1")
    ax.axvline(100, linestyle="--", linewidth=1)
    ax.legend()
    ax.grid(axis="x", alpha=.25)
    _save(fig, out, "Figure5_SIA_reach_requirements")


def figure5b(c, out):
    fig, ax1 = plt.subplots(figsize=(11, 6.5))
    ax2 = ax1.twinx()
    ax1.plot(c.SIA_reach_pct, c.U5_weighted_Rcov, marker="o", label="U5-weighted Rcov")
    ax2.plot(c.SIA_reach_pct, c.provinces_Rcov_lt1, marker="s", linestyle="--", label="Provinces with Rcov < 1")
    ax1.axhline(1, linestyle="--", linewidth=1)
    ax1.set_xlabel("Random SIA reach (% of modeled U5 population)")
    ax1.set_ylabel("U5-weighted Rcov")
    ax2.set_ylabel("Number of provinces with Rcov < 1")
    ax1.set_title("Coverage-adjusted Transmission Potential Under Increasing Random SIA Reach")
    lines = ax1.get_lines()[:1] + ax2.get_lines()[:1]
    ax1.legend(lines, [x.get_label() for x in lines], loc="upper right")
    ax1.grid(alpha=.25)
    _save(fig, out, "Figure5B_random_SIA_scenario")


def figure6(s, out):
    d=s.set_index("province").loc[TABLE3_PROVINCES].reset_index()
    y=np.arange(len(d))

    fig, (ax1,ax2)=plt.subplots(1,2,figsize=(14,6.8))
    ax1.barh(y,100*d.outbreak_probability)
    ax1.set_yticks(y,d.province)
    ax1.invert_yaxis()
    ax1.set_xlabel("Establishment probability (%)")
    ax1.set_title("A. Single-importation establishment probability")
    for yy,val in zip(y,100*d.outbreak_probability):
        ax1.text(val+0.5,yy,f"{val:.1f}%",va="center",fontsize=8)
    ax2.scatter(d.r_cov,100*d.outbreak_probability)
    for row in d.itertuples():
        ax2.annotate(row.province,(row.r_cov,100*row.outbreak_probability),xytext=(4,4),textcoords="offset points",fontsize=8)
    ax2.set_xlabel("Rcov")
    ax2.set_ylabel("Establishment probability (%)")
    ax2.set_title("B. Rcov versus simulated establishment")
    for ax in (ax1,ax2): ax.grid(alpha=.2)
    _save(fig,out,"Figure6A_stochastic_establishment")

    fig,ax=plt.subplots(figsize=(11,6.5))
    x=np.arange(len(d)); width=.38
    ax.bar(x-width/2,d.deterministic_peak,width,label="Deterministic peak I")
    ax.bar(x+width/2,d.peak_median_established,width,label="Stochastic median peak I")
    ax.errorbar(x+width/2,d.peak_median_established,
                yerr=np.vstack([d.peak_median_established-d.peak_p05_established,
                                d.peak_p95_established-d.peak_median_established]),
                fmt="none",capsize=3)
    ax.set_xticks(x,d.province,rotation=35,ha="right")
    ax.set_ylabel("Peak infectious count")
    ax.set_title("Deterministic and Stochastic Peak Infectious Counts")
    ax.legend(); ax.grid(axis="y",alpha=.2)
    _save(fig,out,"Figure6B_deterministic_vs_stochastic")

    fig,ax=plt.subplots(figsize=(10,6.5))
    ax.bar(d.province,100*d.non_establishment_probability)
    ax.set_ylabel("Non-establishment probability (%)")
    ax.set_title("Stochastic Non-establishment Probability Following One Imported Case")
    ax.tick_params(axis="x",rotation=35); ax.grid(axis="y",alpha=.2)
    _save(fig,out,"Figure6C_stochastic_nonestablishment")


def figure7(k,out):
    vals=[
        ("Required immune fraction", P.hit*100),
        ("Theoretical one-dose requirement", float(k["one_dose_required_pct"])),
        ("Complete two-dose-series requirement", float(k["two_dose_required_pct"])),
        ("Current MCV1 (12–23 months)", float(k["direct_national_mics_mcv1_pct"])),
        ("Current MCV2 (24–35 months)", 36.8),
    ]
    labels=[v[0] for v in vals]; y=[v[1] for v in vals]
    fig,ax=plt.subplots(figsize=(10,6))
    ax.barh(labels,y)
    ax.set_xlim(0,105); ax.set_xlabel("Coverage / immune fraction (%)")
    ax.set_title("Theoretical One-dose and Complete Two-dose-series Requirements")
    for i,v in enumerate(y): ax.text(v+0.6,i,f"{v:.1f}%",va="center")
    ax.grid(axis="x",alpha=.2)
    _save(fig,out,"Figure7_single_vs_two_dose")


def _norm_name(x: str) -> str:
    x=unicodedata.normalize("NFKD",str(x)).encode("ascii","ignore").decode().strip()
    aliases={"Jawzjan":"Jowzjan","Day Kundi":"Daykundi","Daikundi":"Daykundi","Sar-e-Pul":"Sar-e Pul","Sar-e Pol":"Sar-e Pul"}
    return aliases.get(x,x)


def maps(p, who, gadm: Path, out: Path):
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise SystemExit("Map regeneration requires geopandas. Install requirements-optional.txt") from exc
    g=gpd.read_file(gadm)
    namecol=next((c for c in ["NAME_1","VARNAME_1","NL_NAME_1"] if c in g.columns),None)
    if not namecol:
        raise ValueError("Could not find a GADM level-1 province-name field (expected NAME_1).")
    g["province"] = g[namecol].map(_norm_name)
    p2=p.copy(); p2["province"]=p2.province.map(_norm_name)
    m=g.merge(p2,on="province",how="left",validate="one_to_one")
    if m.r_cov.isna().any():
        raise ValueError("Unmatched GADM province names: "+", ".join(m.loc[m.r_cov.isna(),"province"]))

    fig,ax=plt.subplots(figsize=(10,8)); m.plot(column="r_cov",legend=True,edgecolor="black",linewidth=.4,ax=ax)
    ax.set_title("Coverage-adjusted Transmission Potential (Rcov) by Province")
    ax.set_axis_off(); _save(fig,out,"Figure7A_Rcov_map")

    w=who.copy(); w["province"]=w.province.map(_norm_name)
    m2=m.merge(w[["province","exact_value_if_stated"]],on="province",how="left")
    fig,(a,b)=plt.subplots(1,2,figsize=(15,7))
    m2.plot(column="r_cov",legend=True,edgecolor="black",linewidth=.35,ax=a); a.set_title("A. Rcov")
    m2.plot(color="lightgrey",edgecolor="black",linewidth=.35,ax=b)
    exact=m2[m2.exact_value_if_stated.notna()]
    exact.plot(column="exact_value_if_stated",legend=True,edgecolor="black",linewidth=.35,ax=b)
    b.set_title("B. Four highest exact WHO week-14 incidence values")
    for ax in (a,b): ax.set_axis_off()
    _save(fig,out,"Figure7B_modelled_vs_observed")

    def cat(v):
        if v>=.8:return "≥80%"
        if v>=.6:return "60–<80%"
        if v>=.4:return "40–<60%"
        return "<40%"
    m2["reach_category"]=m2.sia_preferential.map(cat)
    fig,ax=plt.subplots(figsize=(10,8)); m2.plot(column="reach_category",categorical=True,legend=True,edgecolor="black",linewidth=.4,ax=ax)
    ax.set_title("Idealized Preferential-targeting SIA Reach Category")
    ax.set_axis_off(); _save(fig,out,"Figure7C_preferential_SIA_map")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--tables-dir",type=Path,default=ROOT/"outputs"/"tables")
    p.add_argument("--output-dir",type=Path,default=ROOT/"outputs"/"figures")
    p.add_argument("--gadm",type=Path,default=None,help="Optional GADM v4.1 Afghanistan level-1 boundary file (.gpkg/.shp/.json).")
    args=p.parse_args()
    prov,key,stoch,campaign,traj=_load(args.tables_dir)
    figure1(prov,key,args.output_dir); figure2(prov,key,args.output_dir)
    figure3(traj,args.output_dir); figure4(prov,args.output_dir)
    figure5(prov,args.output_dir); figure5b(campaign,args.output_dir)
    figure6(stoch,args.output_dir); figure7(key,args.output_dir)
    if args.gadm:
        who=pd.read_csv(ROOT/"data"/"derived"/"who_week14_measles_incidence_categories.csv")
        maps(prov,who,args.gadm,args.output_dir)
    else:
        print("Map figures 7A–7C skipped: pass --gadm PATH to a GADM v4.1 Afghanistan level-1 boundary file.")

if __name__=="__main__":
    main()
