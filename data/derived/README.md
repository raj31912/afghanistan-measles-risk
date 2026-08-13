# Derived analysis inputs

These four CSV files are the compact, non-identifiable inputs required by the audited manuscript pipeline. They contain no household-level records.

- `mics_official_sampling_errors.csv`: 34 published provincial MCV1 point estimates, standard errors, design effects, sample denominators and 95% CIs from Afghanistan MICS 2022–23 sampling-error tables. **Row order is preserved from the audited analysis because it controls fixed-seed Monte Carlo draw ordering.**
- `mics_national_official.csv`: direct national design-based MCV1 estimate and uncertainty/sample denominators.
- `population_2026_provincial_clean.csv`: province-level total and under-five population denominators used in the model.
- `who_week14_measles_incidence_categories.csv`: published WHO week-14 ordinal incidence category (1–5) for all 34 provinces plus exact values only where explicitly stated.

See `../../docs/DATA_DICTIONARY.md` and the root `README.md` for provenance and field definitions.
