#!/usr/bin/env python3
"""Fail-fast numerical validation of generated outputs against locked manuscript values."""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
TABLES=ROOT/"outputs"/"tables"

EXPECTED={
"direct_national_mics_mcv1_pct":51.2,
"u5_weighted_r_cov":7.945929252533159,
"routine_mcv1_unprotected_3p5yr_n":3020631.692617184,
"sia_preferential_n":4056042.0399888135,
"sia_random_n":7530647.856211407,
"one_dose_required_pct":100.35842293906809,
"two_dose_required_pct":96.21993127147766,
"week14_kendall_tau_b":0.33104795088367867,
"week14_spearman_rho":0.4302580494917712,
"monte_carlo_r_cov_ci_low":7.690345305309093,
"monte_carlo_r_cov_ci_high":8.202344006034181,
}


def main():
    f=TABLES/"key_results_recomputed.csv"
    if not f.exists():
        raise SystemExit("Run scripts/run_analysis.py first.")
    k=pd.read_csv(f).set_index("metric")["value"]
    failed=[]
    for name,exp in EXPECTED.items():
        got=float(k[name])
        if not np.isclose(got,exp,rtol=0,atol=1e-10): failed.append((name,got,exp))
    s=pd.read_csv(TABLES/"stochastic_results_corrected.csv")
    rng=(100*s.outbreak_probability.min(),100*s.outbreak_probability.max())
    if not (np.isclose(rng[0],65.9) and np.isclose(rng[1],93.9)):
        failed.append(("stochastic establishment range",rng,(65.9,93.9)))
    if failed:
        for x in failed: print("FAIL",x,file=sys.stderr)
        raise SystemExit(1)
    print("PASS: locked headline, Monte Carlo, and stochastic results reproduced.")

if __name__=="__main__": main()
