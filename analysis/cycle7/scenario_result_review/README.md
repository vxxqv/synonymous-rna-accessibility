# Independent finite-scenario verification

Requirements: Python3.12 or newer, standard library only. No installation, Git metadata or network is needed at runtime.

Suggested repository integration path:

`analysis/cycle7/scenario_result_review/verify_observed.py`

Run from the computational repository root with a fresh output directory outside it:

```
python analysis/cycle7/scenario_result_review/verify_observed.py --repo . --out ../NEW_SCENARIO_INDEPENDENT_REVIEW
```

The verifier reads the fixed default files:

- results/cycle6/full_result_v1.json.gz
- results/cycle7/scenario_geometry.json.gz
- analysis/cycle7/freeze.json and its six hash-pinned source files, read only for integrity checks

Alternatively provide --primary PATH, --result PATH and/or --freeze PATH. The two result inputs may be exact plain JSON or gzip; decompressed content hashes are fixed in the verifier, and stored-file hashes are retained in the receipt. Changing serialization bytes or scientific inputs fails. The script rejects existing output, output inside the repository and -O/-OO. It never imports the production algorithms.

Outputs are verification.json and all_pairs.csv, covering all2556 pairs and all69012 capacities. Historical Git access order is not reauthenticated, and historical_git_verified is always false. The separate source_commit_check.json is a review-time content comparison, not a dependency of replay.

For small synthetic controls, copy test_synthetic.py beside verify_observed.py in a fresh directory and run it before creating its synthetic_receipt.json. Those tests contain no biological data. The original synthetic receipt is retained and never overwritten. portable_checks.py and portable_fixture/ are review-workspace QA scaffolding; neither is needed in the repository or CI.

The output quantifies exact descriptive route existence across27 dependent preparation fields under one common absolute threshold. It provides no confidence calibration, independent biological validation, causal mechanism or whole-synonymous-space connectivity claim.
