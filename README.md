# Afghanistan measles outbreak-risk analysis — reproducibility repository

This repository contains the audited analysis pipeline and non-identifiable derived inputs used for the revised Vaccine manuscript:

**Measles Outbreak Risk in Afghanistan: Provincial Immunity Gaps and Epidemic Potential Following Immunization System Disruption**

The repository was rebuilt for the major revision so that the manuscript's numerical results come from one documented pipeline rather than from mixed analysis versions or hard-coded provincial outputs.

## What this repository reproduces

A full run reproduces the revised manuscript's central analytical results, including:

- direct national MICS MCV1 coverage: **51.2%** (95% CI 49.1–53.2; weighted n=6,383; unweighted n=6,177);
- U5-weighted provincial MCV1 aggregate: **50.57%**;
- U5-weighted coverage-adjusted transmission-potential index (Rcov): **7.95** (conditional Monte Carlo 95% interval 7.69–8.20);
- provincial Rcov range: **2.63–13.80**, with all 34 provinces above Rcov=1 under the central assumptions;
- standardized 3.5-year upper-bound routine-MCV1 protection gap: **3.02 million** children (conditional interval 2.92–3.12 million);
- idealized preferential-targeting SIA requirement: **4.06 million** children (**49.8%** of modeled U5);
- random-targeting SIA requirement: **7.53 million** children (**92.4%** of modeled U5);
- theoretical one-dose coverage requirement: **100.36%** at R0=15 and VE1=0.93;
- theoretical complete two-dose-series requirement: **96.22%** at R0=15 and empirical VE2=0.97;
- corrected deterministic SEIR results, structural sensitivity, two-dose VE sensitivity, and natural-history sensitivity;
- stochastic tau-leaping establishment probabilities of **65.9%–93.9%** across nine illustrative provinces, using 1,000 simulations per province, dt=0.1 day and base random seed 42;
- WHO week-14 retrospective ordinal geographic concordance: Kendall tau-b **0.331** (p=0.011) and Spearman rho **0.430** (p=0.011).

The regression test suite protects these values against accidental analytical drift.

## Interpretation guardrails

These distinctions are essential to the manuscript and to correct reuse of the code:

1. **Rcov is not an empirically estimated province-specific effective reproduction number.** It is a coverage-adjusted transmission-potential index under explicit homogeneous-mixing assumptions.
2. The **3.5-year quantity is an upper-bound routine-MCV1 protection gap**, not the actual susceptible population in early 2025. Natural infection, SIAs, MCV2, maternal antibody, migration, and time-varying coverage are not dynamically subtracted.
3. **Preferential-targeting SIA** is an idealized lower bound. **Random-targeting SIA** assumes reach is random within the modeled U5 population. The equations are different and must not be interchanged.
4. Theoretical SIA requirements are **not capped at 100%**. Values above 100% identify one-pass infeasibility under the stated assumptions.
5. Stochastic P05–P95 intervals quantify **demographic stochasticity conditional on fixed inputs**, not MICS survey sampling uncertainty.
6. The WHO comparison is **retrospective geographic concordance, not external model validation**. Published WHO map categories are treated as ordinal categories; map-bin midpoints are never treated as observed numeric incidence.
7. The direct national MICS estimate (**51.2%**) is distinct from the U5-weighted aggregate of provincial estimates (**50.57%**).

## Repository structure

```text
.
├── data/
│   ├── derived/                 # non-identifiable analysis inputs derived from public sources
│   └── external/                # instructions for source data not redistributed here
├── docs/
│   ├── DATA_DICTIONARY.md
│   ├── MANUSCRIPT_OUTPUT_MAP.md
│   ├── NUMERICAL_AUDIT.md
│   ├── REPRODUCIBILITY.md
│   └── RELEASE_CHECKLIST.md
├── outputs/
│   ├── figures/                 # regenerated analytical figures
│   └── tables/                  # exact numerical and manuscript-display CSV outputs
├── scripts/
│   ├── extract_mics_mcv1.py     # optional independent raw-MICS point-estimate check
│   ├── generate_figures.py
│   ├── reproduce_all.py
│   ├── run_analysis.py
│   └── validate_outputs.py
├── src/measles_analysis/
│   ├── core.py
│   ├── pipeline.py
│   └── tables.py
├── tests/
├── requirements.txt
└── pyproject.toml
```

## Data provenance

The pipeline uses four compact, non-identifiable derived inputs in `data/derived/`.

### 1. Afghanistan MICS 2022–23

Public survey information and report:

- MICS survey portal: `https://mics.unicef.org/surveys`
- Afghanistan MICS 2022–23 Survey Findings Report: `https://www.unicef.org/afghanistan/reports/afghanistan-multiple-indicator-cluster-survey-mics-2022-2023`

`mics_official_sampling_errors.csv` contains the published provincial MCV1 point estimates, standard errors and 95% confidence intervals from the MICS sampling-error tables. `mics_national_official.csv` contains the direct national design-based estimate.

**Raw MICS microdata are not redistributed in this repository.** Researchers with authorized access to `ch2022.sav` can run `scripts/extract_mics_mcv1.py` to independently reproduce the MCV1 point estimates. The model itself uses the published official provincial estimates/SEs/CIs, so the raw SAV file is not required to reproduce the manuscript's model results.

### 2. OCHA/HDX population data

Source: `https://data.humdata.org/dataset/cod-ps-afg`

`population_2026_provincial_clean.csv` contains the provincial total and under-five population denominators used in the audited analysis.

### 3. WHO week-14 2025 surveillance

Source report: `https://www.emro.who.int/images/stories/afghanistan/Outbreak-Situation-Report-Week-14-2025.pdf`

`who_week14_measles_incidence_categories.csv` transcribes the five published ordinal provincial incidence categories from WHO Figure 3. Exact numeric values are retained only for the four provinces for which exact values were explicitly reported. Range midpoints are not used as observations.

### 4. GADM boundaries for map regeneration

Source: `https://gadm.org/`

GADM v4.1 Afghanistan level-1 boundaries are **not redistributed**. They are required only to regenerate Figures 7A–7C. See `data/external/README.md` and the `--gadm` option in `scripts/generate_figures.py`.

## Installation

A clean virtual environment is recommended.

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install -r requirements.txt
```

The package was regression-tested under Python 3.13.5 with the versions pinned in `requirements.txt`.

Optional raw-MICS extraction and map regeneration require additional packages:

```bash
pip install -r requirements-optional.txt
```

## Reproduce the analysis

From the repository root:

```bash
python scripts/run_analysis.py
python scripts/validate_outputs.py
pytest -q
```

The full stochastic run is included by default and typically takes only several seconds on a modern desktop. For a quick deterministic/core run:

```bash
python scripts/run_analysis.py --skip-stochastic
```

To regenerate all non-map analytical figures after a full run:

```bash
python scripts/generate_figures.py
```

To regenerate the map figures as well:

```bash
python scripts/generate_figures.py --gadm /path/to/gadm41_AFG_level1.gpkg
```

A convenience command that runs the complete analysis, numerical validation, tests, and non-map figure regeneration is:

```bash
python scripts/reproduce_all.py
```

## Exact reproducibility settings

Central analytical parameters are defined once in `src/measles_analysis/core.py`:

- R0 = 15 (structural sensitivity 12–18)
- VE1 = 0.93
- empirical VE2 = 0.97 (sensitivity 0.96–0.99)
- latent period = 8 days (sensitivity 6–8)
- infectious-compartment duration = 4 days (sensitivity 3–5)
- deterministic horizon = 365 days; output step = 0.25 day
- stochastic simulations = 1,000 per selected province
- tau-leaping dt = 0.1 day
- base stochastic seed = 42, transformed into stable province-specific streams using SHA-256
- outbreak establishment criterion = peak I > 10
- survey-propagation Monte Carlo draws = 100,000
- Monte Carlo seed = 20260810

The order of rows in `mics_official_sampling_errors.csv` is intentionally preserved because the fixed-seed Monte Carlo draws are generated column-by-column in that order. Analytical point estimates are order-invariant; preserving the audited row order makes the conditional Monte Carlo interval exactly reproducible byte-for-byte.

## Outputs

`outputs/tables/` contains both exact numerical results and manuscript-facing display tables. Important files include:

- `key_results_recomputed.csv`
- `provincial_results_recomputed.csv`
- `deterministic_seir_results.csv`
- `stochastic_results_corrected.csv`
- `structural_sensitivity_results.csv`
- `two_dose_ve_sensitivity.csv`
- `seir_natural_history_sensitivity.csv`
- `campaign_scenario_results.csv`
- `table1_manuscript.csv`
- `table2_manuscript.csv`
- `table3_manuscript.csv`
- `analysis_metadata.json`

See `docs/MANUSCRIPT_OUTPUT_MAP.md` for the mapping from each manuscript result/table/figure to its generated source.

## Reproducibility tests

The test suite checks:

- core formulas, including the corrected random-targeting SIA denominator;
- no silent clipping of theoretical SIA requirements above 100%;
- the final six deterministic Figure 3/Table 2 provinces;
- the nine stochastic provinces and final Table 3 order;
- all headline and Monte Carlo values;
- campaign-scenario values at 70%, 90%, 95% and 100% random SIA reach;
- natural-history sensitivity against final Appendix A Table A3;
- locked stochastic probabilities and peak intervals;
- final manuscript table display ordering and rounding.

## Figures and journal styling

`scripts/generate_figures.py` recreates the analytical content of the figures from generated data. Map figures require an externally supplied GADM boundary file. The journal-submission figure files may contain cosmetic layout refinements (font sizing, label placement, scale-bar placement) applied after analytical generation; such cosmetic refinements do not change plotted values. The CSV outputs are the numerical source of truth.

## Privacy and ethics

The repository contains no names, addresses, household identifiers or individual-level survey records. The bundled analysis inputs are province-level aggregate values from public reports/data sources. Raw MICS microdata are intentionally excluded.

## Repository Status 

This repository contains the reproducibility code, derived non-identifiable analysis inputs, and generated outputs associated with the manuscript Measles Outbreak Risk in Afghanistan: Provincial Immunity Gaps and Epidemic Potential Following Immunization System Disruption.

The individual-level Afghanistan MICS 2022–23 microdata are not redistributed in this repository. Researchers with authorized access to the MICS child dataset can use the provided extraction script to independently reproduce the vaccination-coverage estimates.

A permanent DOI will be added following archival release through Zenodo. 
