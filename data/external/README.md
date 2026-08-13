# External source files not redistributed

The core manuscript analysis does **not** require raw individual-level data. The following external files are intentionally not redistributed.

## Optional Afghanistan MICS raw microdata

`ch2022.sav` may be obtained through the UNICEF MICS data-access process (`https://mics.unicef.org/surveys`). If authorized, it can be checked with:

```bash
python scripts/extract_mics_mcv1.py /path/to/ch2022.sav \
  --official data/derived/mics_official_sampling_errors.csv
```

Expected independent audit results:
- 6,177 unweighted children aged 12–23 months;
- weighted denominator approximately 6,383.374;
- national MCV1 approximately 51.1668% (reported 51.2%);
- all 34 raw provincial point estimates within 0.05 percentage points of the published rounded provincial estimates.

The raw SAV is not required by `scripts/run_analysis.py`, which uses the official published provincial estimates and standard errors.

## Optional GADM v4.1 boundaries

Download an Afghanistan level-1 boundary file from `https://gadm.org/` and pass it to `scripts/generate_figures.py --gadm PATH` to regenerate map panels. GADM files are not bundled with this repository.
