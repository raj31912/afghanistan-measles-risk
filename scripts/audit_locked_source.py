#!/usr/bin/env python3
"""Cross-check all locked numerical source-of-truth metrics against generated outputs."""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"outputs"/"tables"
DATA=ROOT/"data"/"derived"


def close(a,b,tol=1e-8): return np.isclose(float(a),float(b),rtol=0,atol=tol)

def main():
    key=pd.read_csv(OUT/"key_results_recomputed.csv").set_index("metric")["value"]
    prov=pd.read_csv(OUT/"provincial_results_recomputed.csv")
    det=pd.read_csv(OUT/"deterministic_seir_results.csv")
    sto=pd.read_csv(OUT/"stochastic_results_corrected.csv")
    struct=pd.read_csv(OUT/"structural_sensitivity_results.csv")
    nat=pd.read_csv(OUT/"seir_natural_history_sensitivity.csv")
    national=pd.read_csv(DATA/"mics_national_official.csv").iloc[0]

    checks=[]
    def add(label,got,expected,tol=1e-8): checks.append((label,got,expected,close(got,expected,tol)))

    add("National MCV1",100*national.mcv1,51.2)
    add("National MCV1 SE",100*national.mcv1_se,1.0)
    add("National weighted n",national.weighted_n,6383)
    add("National unweighted n",national.unweighted_n,6177)
    add("Provincial MCV1 min",100*prov.mcv1.min(),8.6)
    add("Provincial MCV1 max",100*prov.mcv1.max(),88.7)
    add("Direct national Rcov",key.direct_national_r_cov,7.8576)
    add("U5-weighted provincial MCV1",key.u5_weighted_provincial_mcv1_pct,50.56681539402753)
    add("U5-weighted Rcov",key.u5_weighted_r_cov,7.945929252533159)
    add("Rcov MC low",key.monte_carlo_r_cov_ci_low,7.690345305309093)
    add("Rcov MC high",key.monte_carlo_r_cov_ci_high,8.202344006034181)
    add("Provincial Rcov min",prov.r_cov.min(),2.62635)
    add("Provincial Rcov max",prov.r_cov.max(),13.8003)
    add("Provinces Rcov>1",(prov.r_cov>1).sum(),34)

    age=struct[struct.scenario.str.startswith("Age-eligibility")].iloc[0]
    add("Age proxy U5-weighted MCV1",age.u5_weighted_mcv1_pct,42.98179308492339)
    add("Age proxy U5-weighted Rcov",age.u5_weighted_r_cov,9.004039864653187)
    add("National U5",key.national_u5_population,8146035.52681)
    add("Annual routine gap",prov.annual_routine_mcv1_unprotected.sum(),863037.6264620527)
    add("3.5-y routine gap",key.routine_mcv1_unprotected_3p5yr_n,3020631.692617184)
    add("3.5-y gap MC low",key.monte_carlo_gap_ci_low_n,2923471.88328927)
    add("3.5-y gap MC high",key.monte_carlo_gap_ci_high_n,3118107.3315593363)
    add("Preferential SIA n",key.sia_preferential_n,4056042.0399888135)
    add("Preferential SIA %",key.sia_preferential_pct_u5,49.79160754504057)
    add("Preferential SIA MC low",key.monte_carlo_sia_pref_ci_low_n,3906795.0210519265)
    add("Preferential SIA MC high",key.monte_carlo_sia_pref_ci_high_n,4205774.204262473)
    add("Random SIA n",key.sia_random_n,7530647.856211407)
    add("Random SIA %",key.sia_random_pct_u5,92.44555626386301)
    add("Random SIA MC low",key.monte_carlo_sia_rand_ci_low_n,7412834.064133138)
    add("Random SIA MC high",key.monte_carlo_sia_rand_ci_high_n,7586402.775584033)
    add("One-dose requirement",key.one_dose_required_pct,100.35842293906809)
    add("Two-dose requirement",key.two_dose_required_pct,96.21993127147766)
    add("Deterministic total-U5 AR min",100*det.attack_rate_total_u5.min(),15.894277473166188)
    add("Deterministic total-U5 AR max",100*det.attack_rate_total_u5.max(),92.00114726818192)
    add("Deterministic susceptible AR min",100*det.attack_rate_initial_susceptible.min(),90.78517038076232)
    add("Deterministic susceptible AR max",100*det.attack_rate_initial_susceptible.max(),99.99989846588716)
    add("Stochastic establishment min",100*sto.outbreak_probability.min(),65.9)
    add("Stochastic establishment max",100*sto.outbreak_probability.max(),93.9)
    add("Stochastic sims/province",sto.n_sims.min(),1000)
    add("Latent central",8,8)
    add("Latent sensitivity low",nat.latent_days.min(),6)
    add("Infectious central",4,4)
    add("Infectious sensitivity low",nat.infectious_days.min(),3)
    add("Infectious sensitivity high",nat.infectious_days.max(),5)
    add("Approx generation interval",8+4,12)
    add("Week14 Kendall",key.week14_kendall_tau_b,0.33104795088367867)
    add("Week14 Spearman",key.week14_spearman_rho,0.4302580494917712)

    failed=[x for x in checks if not x[3]]
    if failed:
        for label,got,exp,_ in failed: print(f"FAIL {label}: got {got!r}, expected {exp!r}",file=sys.stderr)
        raise SystemExit(1)
    print(f"PASS: {len(checks)} locked numerical checks reproduced.")

if __name__=="__main__": main()
