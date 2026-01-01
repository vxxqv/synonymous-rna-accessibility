# PERSIST revision 2 handoff

READY for independent review. No real sequence has been folded, no real sequence predictor has been calculated, and no cellular target or deposited AUP value has been decoded by this task. The proposed freeze intentionally authorizes no execution. No repository, manuscript or account action was performed.

## Review point

| Artifact | SHA-256 |
|---|---|
| `protocol.json` | `a3dc686a23acaa12af2f2f44154b8baac819cb6df3a0a6d81328dd3ebdc5fd24` |
| `persist_prediction.py` | `962f953ed17574db96f5ba0cde4827486830244941ea745bb34fb697a5514aaa` |
| `test_persist_prediction.py` | `a261667742e82a83328f39940acb998303b82dc7d2434daf0a8ff26dc250934c` |
| `protocol.md` | `26531e3e02493af6554fb3f693a82dbaf4f6fab160a90e30b8df10708504be27` |
| `metadata_settings.json` | `ebc6a30f569ae12271e74837e180a7f9d98d4e222e29a14046c9ec007e455935` |
| `synthetic_test_receipt.json` | `fad868cfbed9ce423be4f8c260654dfe047fdcc0ee50a8aaaa44947ccc4abaaa` |
| `proposed_freeze.json` | `6ea3028cf4596ad062c536874749d9a4781cfeffd4e3127673303f9836390782` |

The full source and parameter manifest is embedded in the protocol and proposed freeze. All 37 identifiers, 22/15 allocation, 15 development groups and per-input hashes are preserved. The runtime is Python 3.12.14, NumPy 2.5.4 and ViennaRNA 2.7.2. Every real stage checks runtime versions and input hashes. Effective binding-exposed ensemble settings and exported Turner2004 parameters are included. The unexposed pair-dissimilarity field keeps the pinned library default; the runner never calls comparative folding.

## Completed checks

The final synthetic suite passed all 47 tests in 1.26 seconds of pytest time, 1.80 seconds wall time. Source hashes remained unchanged through the run. The command-line help check passed. Tests cover one-based upper-triangle indexing, double counting and both MUP identities; probability and finite-value gates; the installed 12-position synthetic folding interface; ridge agreement with an augmented least-squares oracle; constant and collinear columns; seven positive response scale/translation cases; homogeneous penalty ties; exact zero-range and zero-baseline cases; poisoned held-out target cells; training-only fold preparation; malformed, missing and duplicate targets; feature-artifact identity; prediction hash and committed-blob gates; runtime checks; and complete influence reporting.

Two preflight failures are retained in `preflight_incidents.json`. A metadata assertion initially compared RNA and DNA alphabets literally; applying the already authorized T-to-U normalization resolved it without changing a row. The initial synthetic folding test found an unsupported manual `pf_scale` setter. The final implementation uses the supported `exp_params_rescale(mfe)` API once and records its result. The original settings snapshot is preserved as `metadata_settings_initial.json`. Neither failure involved a real fold or cellular outcome access.

## Execution interface after root authorization

From the workspace root, the following syntax has been checked against the CLI help. These commands have not been run. The execution freeze paths and output locations are examples for root to create after review. Output directories must not exist. Never change a frozen file in place to accommodate an observed failure.

```powershell
& 'work/env/Scripts/python.exe' 'work/persist_cycle3/prediction_v2/persist_prediction.py' features --protocol 'work/persist_cycle3/prediction_v2/protocol.json' --freeze 'work/persist_cycle3/prediction_v2/root_features_freeze.json' --root '.' --output 'work/persist_cycle3/prediction_v2/features_run'
& 'work/env/Scripts/python.exe' 'work/persist_cycle3/prediction_v2/persist_prediction.py' train --protocol 'work/persist_cycle3/prediction_v2/protocol.json' --freeze 'work/persist_cycle3/prediction_v2/root_train_freeze.json' --root '.' --output 'work/persist_cycle3/prediction_v2/training_run' --features 'work/persist_cycle3/prediction_v2/features_run/features.json'
& 'work/env/Scripts/python.exe' 'work/persist_cycle3/prediction_v2/persist_prediction.py' validate --protocol 'work/persist_cycle3/prediction_v2/protocol.json' --freeze 'work/persist_cycle3/prediction_v2/root_validation_freeze.json' --root '.' --output 'work/persist_cycle3/prediction_v2/validation_run' --predictions 'work/persist_cycle3/prediction_v2/training_run/predictions.json' --prediction-freeze 'work/persist_cycle3/prediction_v2/root_prediction_commit.json'
```

An authorized execution freeze needs `status: authorized_by_root`, the relevant `authorized_stages` list, and `input_sha256` containing at least the protocol and module plus all desired pinned inputs. The protocol itself checks its complete source manifest. The training freeze must additionally hash `features.json`. Root should authorize stages separately. The current candidate freeze uses `status: proposed_for_root_review` and an empty stage list, so it fails closed.

The separate prediction receipt requires `prediction_sha256`, `repository`, `commit`, and repository-relative `git_path`. Validation verifies both the local bytes and `git show commit:git_path` before decoding any holdout response. Root is responsible for making and recording the actual freeze and prediction commits. A committed copy may live in the computational repository while the working copy remains under this directory; its bytes must match exactly.

## Scope of any eventual result

The old deposited AUP provenance failure remains unresolved. The new predictor has a complete explicit computational definition but does not measure cellular structural occupancy. The primary question is incremental prediction of the deposited cellular coefficient on the selected design-source transfer, conditional on its processing and unresolved absolute calibration. Selection, one protein/context, source differences, barcode effects, related constructs and unresolved preparation identity remain limitations. No confidence intervals, biological noise ceiling, calibrated half-life, causal structure claim or accessibility-barrier claim is justified by this pilot.

The next action is independent review of this exact package. Real feature calculation and target access remain stopped pending root's freeze and separate execution authorization.
