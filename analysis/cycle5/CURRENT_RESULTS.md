# Current results and portable verification

The original README is preserved byte-for-byte in its pre-outcome freeze. This separate guide records the completed pilot and post-outcome reproducibility additions.

## Completed pilot and selected-export verification

Both arms are now exposed. `results/cycle5/pilot_closure_report.md`, `pilot_decision.json` and `pilot_evidence_manifest.json` retain the complete failure and descriptive evidence. Original before-outcome labels remain historical. Reserved coverage and the continuous prediction screen fail; no successful independent biological validation is claimed.

The frozen independent checker can verify each selected export without opening the raw biological archive. From the checkout root, use a new output directory:

```text
python -B analysis/cycle5/observed_result_review/verify_observed.py --authorization analysis/cycle5/cambray_pilot/authorization_development_independent_check.json --root . --out NEW_DEVELOPMENT_REVIEW
python -B analysis/cycle5/observed_result_review/verify_observed.py --authorization analysis/cycle5/cambray_pilot/authorization_reserved_validation_independent_check.json --root . --out NEW_RESERVED_REVIEW
```

The independent supplement review scripts preserve their original workspace-relative locations and write beside themselves. Their historical direct execution requires the original disposable workspace layout; do not run them in the immutable results tree or describe them as portable commands. The core checker commands above support explicit roots and output locations. Exact computational agreement does not reconstruct author preprocessing, historical access order or calibrated assay uncertainty.

The reserved supplement review has an explicit scope-only erratum in `results/cycle5/supplement_scope_erratum.json`. Its `correction_v2` package supersedes an incorrect copied development-only scope sentence. Every scientific output and numerical helper is unchanged; original records remain available. The relocation wrapper also retains historical workspace paths and is not a portable checkout command.

A separate location-only adapter now provides portable replay of both byte-pinned independent supplement checkers. It changes no numerical function and requires complete result equality, permitting only a Python-version string difference. Use a new output directory outside the checkout:

```text
python -I -B analysis/cycle5/replay_supplement_reviews.py --repo . --arm development --out ../NEW_EXTERNAL_DEV_REVIEW
python -I -B analysis/cycle5/replay_supplement_reviews.py --repo . --arm reserved_validation --out ../NEW_EXTERNAL_RESERVED_REVIEW
```

Both commands reproduced their authenticated reference JSON bytes under CPython 3.12.14. Receipts are in `results/cycle5/portable_review_checks`. The adapter is a post-outcome reproducibility addition, not a changed analysis or new independent experiment.
