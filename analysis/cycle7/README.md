# Exploratory FAS robustness and joint-scenario verification

The original scientific primary is in `analysis/cycle6` and `results/cycle6`.
This follow-up was frozen at revision `7eae1381c781c93e324f01d4e0d07692bceabea5`
after the primary results had been exposed. It is not fresh validation.
No primary threshold, genotype, normalization, model or missingness rule changed.

The complete result files are `results/cycle7/robustness.json.gz` and
`results/cycle7/scenario_geometry.json.gz`. Their decompressed SHA-256 values are
`861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa` and
`0c2ad1a6a602f80e44ff05eaf3321522d9eaf0940382c599309a04d33f4b383e`.
They retain exact certificates, feasible boundary fields, all omissions,
all 27 scenarios and all 2,556 endpoint pairs. The scientific report and
lossless CSV exports are in `results/cycle7`.

## Independent replay

Use Python 3.12 and new output directories outside the checkout:

```text
python analysis/cycle7/robustness_result_review/verify_portable.py --repo /path/to/checkout --out /new/external/robustness-check
python analysis/cycle7/scenario_result_review/verify_observed.py --repo /path/to/checkout --out /new/external/scenario-check
```

These standard-library implementations require neither Git metadata, network,
SciPy nor the production optimizer. They verify input hashes before processing.
The first reconstructs all 26 sharp radii, 72 observation omissions and 220
stored search certificates. The second independently reconstructs 69,012
scenario capacities, 2,556 complete pair records and their exact threshold
categories. Both checks run in continuous integration.

## Reproduction and figures

The frozen implementation can regenerate the derived exploratory quantities:

```text
python analysis/cycle7/robustness.py --core analysis/cycle6/fas_pilot/ordinal_geometry.py --result results/cycle6/full_result_v1.json.gz --freeze analysis/cycle7/freeze.json --out /new/external/robustness.json
python analysis/cycle7/run_scenarios.py --help
python analysis/cycle7/render_figure.py --robustness results/cycle7/robustness.json.gz --geometry results/cycle7/scenario_geometry.json.gz --out /external/figures
```

Regeneration needs the pinned scientific dependencies. The figure renderer
reads only authenticated results and writes exports outside the computational
archive. No manuscript or publishing assets are included here.

The native reference is fixed at one. The multiplicative boundary concerns
one model-compatible completion of an empirical box, not measured bias or
one correction explaining every preparation. Preparation envelopes are not
confidence intervals. Omitting preparation scenarios narrows them and is
an influence diagnostic, not independent validation. Row 46 is a normalization
identity, not an independently noisy observation. Existing published releases
are preserved; this repository update does not itself create a Zenodo release.
