#!/usr/bin/env python3
"""Run the full audited analysis from repository-relative inputs."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from measles_analysis.pipeline import run_pipeline


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",type=Path,default=ROOT/"data"/"derived")
    p.add_argument("--output-dir",type=Path,default=ROOT/"outputs"/"tables")
    p.add_argument("--skip-stochastic",action="store_true",help="Skip the 9-province x 1,000-run stochastic analysis.")
    args=p.parse_args()
    key=run_pipeline(args.data_dir,args.output_dir,run_stochastic=not args.skip_stochastic)
    print("\nHeadline outputs")
    for k in [
        "direct_national_mics_mcv1_pct","u5_weighted_r_cov","routine_mcv1_unprotected_3p5yr_n",
        "sia_preferential_n","sia_random_n","one_dose_required_pct","two_dose_required_pct",
        "week14_kendall_tau_b","week14_spearman_rho",
    ]:
        print(f"{k}: {key[k]}")


if __name__ == "__main__":
    main()
