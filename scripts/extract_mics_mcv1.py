"""Reproduce Afghanistan MICS 2022-23 MCV1 estimates from ch2022.sav.

This script is intentionally separate from the transmission model so the survey
extraction can be audited independently. It reproduces the published national
MCV1 estimate (51.2%) and the published provincial point estimates.

Required package: pyreadstat
Input: ch2022.sav from the Afghanistan MICS 2022-23 microdata release.

MCV1 construction used here (children age 12-23 months):
  * evidence of first measles vaccination in IM6M1D (valid non-missing value,
    excluding 0 = not given and 99 = no response), OR
  * caretaker report that the child ever received measles vaccination (IM26=1).
The analysis weight is CHWEIGHT and province is HH7.

This definition was checked against the official MICS TC.10/TC.1 vaccination
results. Do not replace it with the older restrictive rule that required IM5=1
and a valid IM6M1Y value; that rule produced materially incorrect provincial
coverage estimates.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def weighted_mean(x: pd.Series, w: pd.Series) -> float:
    mask = x.notna() & w.notna() & (w > 0)
    if not mask.any():
        return np.nan
    return float(np.average(x.loc[mask].astype(float), weights=w.loc[mask].astype(float)))


def normalise_province(name: str) -> str:
    mapping = {
        "Maidan Wardak": "Wardak", "Paktya": "Paktia", "Panjsher": "Panjshir",
        "Jawzjan": "Jowzjan", "Urozgan": "Uruzgan", "Nooristan": "Nuristan",
        "Kunarha": "Kunar", "Sar-e-Pul": "Sar-e Pul",
    }
    return mapping.get(name, name)


def extract_mcv1(sav_path: Path) -> tuple[pd.DataFrame, dict]:
    try:
        import pyreadstat
    except ImportError as exc:
        raise SystemExit(
            "pyreadstat is required. Install with: pip install pyreadstat"
        ) from exc

    usecols = ["HH7", "CAGE", "CHWEIGHT", "IM6M1D", "IM26"]
    df, meta = pyreadstat.read_sav(str(sav_path), usecols=usecols, apply_value_formats=False)

    domain = df.loc[df["CAGE"].between(12, 23, inclusive="both")].copy()
    if len(domain) != 6177:
        raise ValueError(f"Expected 6,177 unweighted children age 12-23 months; found {len(domain):,}.")

    card_or_document = domain["IM6M1D"].notna() & ~domain["IM6M1D"].isin([0, 99])
    caretaker_report = domain["IM26"].eq(1)
    domain["mcv1"] = (card_or_document | caretaker_report).astype(int)

    province_labels = meta.variable_value_labels.get("HH7", {})
    if not province_labels:
        raise ValueError("Province value labels for HH7 were not available in the SAV metadata.")
    domain["province"] = domain["HH7"].map(province_labels).astype(str).map(normalise_province)

    rows = []
    for province, g in domain.groupby("province", sort=True):
        rows.append({
            "province": province,
            "mcv1": weighted_mean(g["mcv1"], g["CHWEIGHT"]),
            "unweighted_n": int(len(g)),
            "weighted_n": float(g["CHWEIGHT"].sum()),
        })
    out = pd.DataFrame(rows).sort_values("province").reset_index(drop=True)
    out["mcv1_pct"] = 100 * out["mcv1"]

    national = {
        "mcv1": weighted_mean(domain["mcv1"], domain["CHWEIGHT"]),
        "mcv1_pct": 100 * weighted_mean(domain["mcv1"], domain["CHWEIGHT"]),
        "unweighted_n": int(len(domain)),
        "weighted_n": float(domain["CHWEIGHT"].sum()),
    }
    return out, national


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("sav", type=Path, help="Path to ch2022.sav")
    p.add_argument("--official", type=Path, default=None,
                   help="Optional mics_official_sampling_errors.csv for cross-checking")
    p.add_argument("--out", type=Path, default=Path("mics_mcv1_from_microdata.csv"))
    args = p.parse_args()

    province, national = extract_mcv1(args.sav)
    province.to_csv(args.out, index=False)
    print(f"National MCV1: {national['mcv1_pct']:.2f}%")
    print(f"Unweighted n: {national['unweighted_n']:,}")
    print(f"Weighted denominator: {national['weighted_n']:.3f}")

    if args.official:
        official = pd.read_csv(args.official)[["province", "mcv1"]]
        check = province.merge(official, on="province", suffixes=("_raw", "_official"), validate="one_to_one")
        check["abs_diff_pp"] = 100 * (check["mcv1_raw"] - check["mcv1_official"]).abs()
        worst = check["abs_diff_pp"].max()
        print(f"Largest absolute difference from official province point estimate: {worst:.3f} percentage points")
        if worst > 0.15:
            raise ValueError("Raw-data extraction does not reproduce official MICS provincial estimates closely enough.")


if __name__ == "__main__":
    main()
