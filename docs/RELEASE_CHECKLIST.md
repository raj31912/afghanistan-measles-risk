# Public-release checklist

Analytical package status: **scientifically reproducible and regression-tested**.

Before a public GitHub/Zenodo release, complete these administrative items:

- [ ] Confirm final repository title and author/contributor list.
- [ ] Choose a code license appropriate for the authors and compatible with source-data terms.
- [ ] Create the public repository/Zenodo record.
- [ ] Record the permanent URL and/or DOI.
- [ ] Run `python scripts/reproduce_all.py` in a clean environment.
- [ ] Confirm `pytest -q` passes.
- [ ] Confirm `python scripts/validate_outputs.py` reports PASS.
- [ ] Upload the repository ZIP or connect Zenodo to the GitHub release.
- [ ] Replace the provisional manuscript Data Availability wording with the actual permanent record.
- [ ] Do not claim a DOI before Zenodo has minted one.

The raw MICS microdata and GADM boundary files should remain excluded from the repository; access instructions are provided instead.
