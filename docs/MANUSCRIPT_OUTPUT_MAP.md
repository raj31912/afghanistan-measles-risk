# Manuscript-to-code output map

This document identifies the numerical source for each major revised-manuscript result.

| Manuscript item | Code/output source |
|---|---|
| National MCV1 51.2% and survey n | `data/derived/mics_national_official.csv`; optional independent check `scripts/extract_mics_mcv1.py` |
| Provincial MCV1 and MICS uncertainty | `data/derived/mics_official_sampling_errors.csv` |
| Rcov formula and provincial Rcov | `src/measles_analysis/core.py`; `outputs/tables/provincial_results_recomputed.csv` |
| U5-weighted Rcov 7.95 and MC interval | `outputs/tables/key_results_recomputed.csv` |
| Population-immunity threshold, one-dose and complete two-dose requirements | `src/measles_analysis/core.py`; `outputs/tables/key_results_recomputed.csv` |
| Deterministic SEIR results | `src/measles_analysis/core.py`; `outputs/tables/deterministic_seir_results.csv` |
| Figure 3 six-province trajectories | `outputs/tables/figure3_seir_trajectories.csv`; `scripts/generate_figures.py` |
| 3.5-year upper-bound routine-MCV1 protection gap | `outputs/tables/provincial_results_recomputed.csv`; `key_results_recomputed.csv` |
| Preferential and random SIA requirements | `core.py`; `provincial_results_recomputed.csv`; `key_results_recomputed.csv` |
| Figure 5B random-reach scenario | `outputs/tables/campaign_scenario_results.csv` |
| Stochastic establishment/non-establishment and peak intervals | `outputs/tables/stochastic_results_corrected.csv` |
| Structural sensitivity / Appendix A Table A1 | `outputs/tables/structural_sensitivity_results.csv` |
| Empirical VE2 sensitivity / Appendix A Table A2 | `outputs/tables/two_dose_ve_sensitivity.csv` |
| Natural-history sensitivity / Appendix A Table A3 | `outputs/tables/seir_natural_history_sensitivity.csv` |
| WHO week-14 Kendall/Spearman concordance | `data/derived/who_week14_measles_incidence_categories.csv`; `key_results_recomputed.csv` |
| Final Table 1 | `outputs/tables/table1_manuscript.csv` |
| Final Table 2 | `outputs/tables/table2_manuscript.csv` |
| Final Table 3 | `outputs/tables/table3_manuscript.csv` |
| Figures 1–7 analytical regeneration | `scripts/generate_figures.py` |
| Figures 7A–7C maps | `scripts/generate_figures.py --gadm ...` plus externally supplied GADM v4.1 level-1 boundary file |

## Locked province sets

Deterministic Figure 3/Table 2:

`Uruzgan, Paktika, Helmand, Kabul, Herat, Bamyan`

Stochastic analysis:

`Uruzgan, Paktika, Nuristan, Herat, Helmand, Kabul, Panjshir, Nimroz, Bamyan`

Final Table 3 display order:

`Uruzgan, Paktika, Nuristan, Helmand, Nimroz, Panjshir, Kabul, Herat, Bamyan`

The distinction between stochastic execution order and Table 3 display order is intentional and regression-tested.
