# Independent selected-result verifier, pre-outcome freeze

Status: synthetic-tested and ready for the primary agent's freeze. No observed development or validation result has been opened. This folder does not authorize observed access.

`verify_observed.py` imports only the hash-pinned independent `reference.py` and Python's standard library. It never imports the primary runner or graph rules. The separate synthetic test does execute the held primary runner to generate invented comparison exports. Neither program discovers biological output paths. The verifier has no archive, network or Git interface.

## Invocation and authorization

```text
work/env/Scripts/python.exe work/joint_function_cycle5/observed_result_review/verify_observed.py --authorization EXPLICIT_AUTHORIZATION.json --root EXPLICIT_PATH_ROOT --out NEW_OUTPUT_DIRECTORY
```

The primary agent must supply a JSON object with schema `Cambray_independent_result_v1`, `execution_authorized: true`, arm exactly `development` or `reserved_validation`, and data_scope exactly `synthetic` or `explicitly_authorized_observed`. Its inputs object must have exactly ten roles: evaluation, unsealed_inputs, predictions, prediction_seal, known_inputs, split, components, edges, adapter, reference. Each role has explicit path and sha256 strings. Relative input paths resolve against the declared root. There are no default input paths, and extra roles, including an archive role, are refused. The adapter role must identify the executing source, and reference must match SHA256 `64a5e9040a77d81dc779583d1d3ea5894a5aebd5718dbf3466281803f48ce375`.

Every input hash is verified before the comparison parses any selected outcome export. The seal is then checked against the declared arm, frozen base hashes and both numerical payloads. Known-input rows and unsealed rows must have their independently derived ID/column whitelists. The seal's arm check precedes decoding the selected export JSONs; each decoded export's arm must also match. This relies on root authorization of the exact file hashes and is not a replacement for the original selective source reader. Input hashes are rechecked at completion. Input files must remain held unchanged to avoid file-system races.

The existing output directory is refused before authorization access. `-O` and `-OO` are refused before imports or input/output. A comparison mismatch stops immediately and writes failure.json rather than a pass receipt. A source-mean hold is verified as `verified_source_hold` without deriving scientific summaries. A successful complete comparison writes receipt.json as `verified_complete_arm`, retaining independent metrics and group summaries. None of these statuses changes the source result's scientific pass/fail interpretation.

## Bounded coverage

The adapter reconstructs all required means, source-mean tolerance decisions and missingness reasons; all original descriptive and prediction exclusions; all three conditional predictions and their withheld-intermediate classification; individual case errors, equal-group MAE/MSE/RMSE and misclassification; pooled confusion/balanced accuracy; every predictive hurdle and failure reason; every complete square's two primary diagonals, raw/log additive projection and native shortfall; point, twelve whole-series, range and numerical-tolerance decisions on both the square and its observed component; and all group/seed/component-stratum primary counts and milestones.

Unknown traits remain None, and known failure in another trait blocks optimistic joint membership. The empirical range threshold is computed from endpoint lower bounds independently. Point/preparation/range numerical tolerances use the relevant full observed component. Thresholds, tolerance values, statuses, labels and coverage decisions compare exactly. Only arithmetic identities such as errors, contrasts and aggregate means have the prespecified absolute/relative comparison tolerance2e-12. That tolerance never changes an allowed vertex, a scientific threshold or a required success condition.

All exact critical geometry cells are also checked: rational boundaries, atoms, intervening open cells, the unbounded first cell, endpoint caps and every represented classification. The frozen reference accepts only int/float point inputs. Consequently the adapter uses Fraction comparisons to construct exact allowed sets, then invokes the same independent transitive closure for geometry. No reference source or scientific rule is changed.

Coverage excludes supplemental amendment diagnostics such as leave-one-group-out and target-balanced recaps; those require their separate verifier. It does not authenticate raw assay processing, historical access order or biological calibration. The reader's stored claims about unparsed source fields are checked for consistency, not independently proven by this post-export replay.

## Synthetic evidence and preserved development history

The final15 tests pass on invented12-square/eight-group fixtures with all four target orientations, one partially measured component-only vertex, exact geometry, shared-edge squares containing multiple globally masked corners, missing known/target values, empty metrics, one-class/zero-error baselines and source provenance holds. Tampered point labels, geometry, per-group errors and sealed payloads fail. Hash, arm, duplicate-key, undeclared-role, overwrite and optimization guards fail as intended. No actual source or outcome input is accepted by the test harness.

```text
work/env/Scripts/python.exe work/joint_function_cycle5/observed_result_review/test_adapter.py
```

This command refuses to overwrite its retained synthetic receipt. For a fresh rerun, copy the owned test/code files into a new review location with the same explicitly required frozen fixture/reference layout, or have the primary agent select a new receipt path without modifying scientific code.

The first synthetic run found an adapter type-bridge bug: Fraction geometry was sent through the frozen reference's intentionally int/float-only finite guard and incorrectly returned an ineligible-endpoint status. The failed receipt is preserved as `synthetic_receipt_v1_type_bridge_failure.json`. Only the adapter's exact geometry bridge was corrected. The subsequent basic and overlap passes are preserved separately. No observed data or mathematical result was used to tune the adapter. `freeze.json` pins the final source, tests, reference and synthetic history before any observed comparison.
