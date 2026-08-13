# Reproducibility notes

## Single-source pipeline

`src/measles_analysis/pipeline.py` is the numerical source of truth for downstream model outputs. It validates the four derived inputs, merges all 34 provinces one-to-one, computes deterministic/core results, propagates sampling uncertainty, runs sensitivity analyses, optionally runs stochastic simulations, and writes machine-readable outputs.

## Survey extraction versus model inputs

The optional raw-MICS extraction script independently checks the MCV1 **point estimates** from authorized `ch2022.sav` microdata. The main model uses official published point estimates, standard errors and confidence intervals from the MICS report. This is deliberate: the public aggregate tables are sufficient for the manuscript model and avoid redistributing individual-level survey records.

## Conditional Monte Carlo sampling uncertainty

Province-summed uncertainty is propagated with 100,000 independent truncated-normal draws using each province's official MCV1 point estimate and SE. Draws are restricted to [0,1]. This is a transparent conditional approximation, not a full multivariate design-based variance estimator.

Because fixed-seed random draws are generated sequentially by province, the audited order of `mics_official_sampling_errors.csv` is preserved. Reordering rows leaves analytical point estimates unchanged but changes the exact fixed-seed Monte Carlo realization and therefore can alter the final digits of the conditional interval.

## Deterministic SEIR

The model uses a U5 population denominator, R0=15, beta=R0×gamma, an 8-day mean latent period and a 4-day infectious-compartment duration. One imported infectious case is removed from the initially susceptible pool, preserving total population mass. Attack rate is computed from susceptible-compartment depletion after introduction.

The deterministic model is an illustrative within-U5 scenario; it is not a province-wide epidemic forecast.

## Stochastic tau-leaping

Nine prespecified illustrative provinces are run with 1,000 simulations each, dt=0.1 day and a 365-day horizon. A stable SHA-256-derived province-specific RNG seed is generated from base seed 42. This avoids accidentally replaying the same random stream in every province while preserving exact reproducibility.

Establishment is operationally defined as peak I>10. Peak P05–P95 intervals are calculated among established outbreaks only. Compartment mass balance is checked during every run.

## Natural-history sensitivity correction

The deposit pipeline uses the same six provinces as final Table 2/Figure 3 and final Appendix A Table A3:

`Uruzgan, Paktika, Helmand, Kabul, Herat, Bamyan`.

An earlier development script still contained an obsolete natural-history sensitivity list including Ghazni and Balkh. That stale list is not used here. A regression test explicitly protects the final six-province specification and Appendix A ranges.

## Table display regression

The pipeline generates manuscript-facing Table 1–3 CSV files with the exact final ordering and rounding. These are regression-tested separately from Word formatting so a cosmetic document edit cannot silently change numerical presentation.

## Figure regeneration

The figure script uses the generated CSV outputs rather than recalculating model results independently. This prevents plots and tables from diverging. GADM boundaries are an external cartographic dependency and are intentionally not redistributed.
