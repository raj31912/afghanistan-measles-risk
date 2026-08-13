# Data dictionary

## `data/derived/mics_official_sampling_errors.csv`

| Field | Meaning |
|---|---|
| `province` | Province name normalized to the 34-province analysis nomenclature. |
| `mcv1` | Official provincial MCV1 point estimate as a proportion (0–1). |
| `mcv1_se` | Official design-based standard error as a proportion. |
| `design_effect` | Published design effect for the MCV1 estimate. |
| `weighted_n` | Published weighted denominator. |
| `unweighted_n` | Published unweighted denominator. |
| `mcv1_ci_low`, `mcv1_ci_high` | Published 95% confidence interval limits as proportions. |
| `mcv1_pct` | MCV1 point estimate in percent; convenience/display field. |
| `source` | Source/provenance note. |

## `data/derived/mics_national_official.csv`

The same core MCV1 fields for the direct national design-based survey estimate. This value is **not** replaced by a U5-population-weighted average of provincial estimates.

## `data/derived/population_2026_provincial_clean.csv`

| Field | Meaning |
|---|---|
| `province` | Province name. |
| `population_total_2026` | Provincial total-population estimate used in Table 1. |
| `population_u5` | Modeled under-five population (0–59 months); SEIR denominator and campaign-scaling denominator. |
| `source` | Source/provenance note. |

## `data/derived/who_week14_measles_incidence_categories.csv`

| Field | Meaning |
|---|---|
| `province` | Province name. |
| `incidence_category` | Ordinal category 1–5 transcribed from WHO week-14 Figure 3. Used for Kendall/Spearman rank comparison. |
| `published_incidence_range_per_10000` | Published legend range corresponding to the category. Used for documentation only. |
| `exact_value_if_stated` | Exact cumulative incidence only where WHO explicitly reported the value; blank otherwise. |
| `source` | WHO report/provenance note. |

**Important:** range midpoints are never converted into numeric observations.

## Main generated fields in `provincial_results_recomputed.csv`

| Field | Meaning |
|---|---|
| `immune_fraction_mcv1` | MCV1 × VE1. |
| `r_cov` | R0 × (1 − MCV1 × VE1); coverage-adjusted transmission-potential index. |
| `sia_preferential` | Idealized preferential-targeting SIA reach fraction. |
| `sia_random` | Random-U5-targeting SIA reach fraction. |
| `*_unattainable_one_pass` | Boolean flag where theoretical reach requirement exceeds 100%. |
| `annual_routine_mcv1_unprotected` | Standardized annual upper-bound routine-MCV1 protection gap. |
| `routine_mcv1_unprotected_3p5yr` | The annual quantity × 3.5 years. |
| `r_cov_ci_low/high` | Province-level propagation of official MICS 95% coverage CI through the monotonic Rcov function. |
| `sia_pref_ci_low/high` | Province-level propagation through preferential SIA formula. |
| `sia_rand_ci_low/high` | Province-level propagation through random SIA formula. |
| `routine_gap_ci_low/high` | Province-level propagation through the upper-bound routine-gap formula. |
| `peak_infectious` | Deterministic maximum concurrent infectious compartment. |
| `peak_day` | Deterministic day of maximum I. |
| `incident_infections_after_introduction` | S(0) − S(T); imported seed excluded. |
| `attack_rate_total_u5` | Incident infections / total modeled U5 population. |
| `attack_rate_initial_susceptible` | Incident infections / initially susceptible population after seed placement. |
