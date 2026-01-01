# Explicit exploratory-successor verifier copy

Synthetic-only readiness copy. No actual run or actual-output authorization is included. The original delivered verifier/ZIP is preserved unchanged. This copy is valid only for the reviewed exploratory operational-cohort successor, not a repaired or successful version1 reproduction.

Only five literals in `verify.py` differ from the original delivery: author/protocol pins, phase VERSION, actual prediction authorization status and operational commit-receipt status. See `adaptation_manifest.json` and `explicit_verifier_changes.diff`. All numerical reference functions, core bytes, numerical tolerances, fixture generator and guard fixtures are unchanged.

- Amendment SHA256: `219dda125d7a8167b2a45fa67499ac4dc24118e6edcbdee83b540f8b77e27f22`.
- Author predict.py SHA256: `be75c90e4076f29952553fa2a22529c2640db0691cc8237e5078c1d1e2cc25d7`.
- This verify.py SHA256: `896f580c33b8f24e3f44180501c51fd141c2b379b1d6f0631aecf5a75558c246`.
- Independent core SHA256: `a5255167a8182d36e107db804d42b87a488e3dc22361cc75c57772a7729c8786`.

Keep `verify.py` adjacent to `reference_core.py`. Use installed CPython3.12.14, NumPy2.5.4, Biopython1.88. There is no network, Git, workbook reader, workspace path or runtime installation. Author source is hashed but never imported by the verifier. The optional synthetic fixture generator imports it only to produce SYN comparison artifacts.

```text
python -I -B verify.py --training-features TRAIN.json --training-targets TRAIN_TARGETS.json --holdout-features HOLDOUT.json --design-receipt DESIGN.json --predictions PREDICTIONS.json --protocol AMENDMENT.md --author-code PREDICT.py --authorization VERIFY_AUTH.json --out NEW_DIRECTORY
```

Add all three scoring options together: `--evaluation EVALUATION.json --holdout-targets HOLDOUT_TARGETS.json --commit-receipt COMMIT_RECEIPT.json`.

The separately frozen verifier authorization must identify `synthetic_only` or, only after separate root permission, `root_authorized_actual_output_verification`; pin the exact verifier/core hashes; and supply an exact `inputs` mapping of every supplied file role to SHA256. The real author artifact must be `root_authorized_operational_prediction`. Its source design receipt must preserve the reviewed operational gate and failed original-reproduction semantics. The evaluator receipt must be `root_verified_operational_prediction_commit` with verified commit and exact prediction/heldout-target hashes. The verifier does not reauthenticate historical Git or recompute source eligibility/protein grouping. Root must authenticate all upstream evidence. Fresh output paths are mandatory; failures remain retained.

Coverage remains9 predictions,35 OOF vectors,7 penalty selections, complete metrics/influence/rate arithmetic, frozen row/error outputs and provenance. Prediction tolerance is1e-9 absolute plus1e-10 relative; score tolerance1e-12 absolute plus1e-10 relative. Discrete decisions use canonical frozen prediction values, and kinetic sign counts describe canonical archived row values after independent numerical arithmetic checks. No scientific5% threshold or estimator changes.

Synthetic replay:2642 checks,13702 numerical values, maximum difference3.921307722976053e-13. All model and scoring payloads are exactly equal to the original algorithm on identical320-row synthetic inputs. Six refusals and37 reciprocal-zero checks passed.

To reproduce readiness in fresh directories, using the exact archived candidate source/amendment:

```text
python -I -B make_synthetic_fixture.py --author-code PREDICT.py --protocol AMENDMENT.md --out NEW_SYNTHETIC_DIRECTORY
python -I -B test_refusals.py --fixture NEW_SYNTHETIC_DIRECTORY --out NEW_REFUSAL_DIRECTORY
```

This verifies a declared calculation only. Any actual result remains post-exposure exploratory prediction within the operational/reference-eligible population, with unavailable true isoforms, decay uncertainty, growth/fit controls and biological independence. No mechanistic or causal calibration, independent experiment or original-gate success is implied.
