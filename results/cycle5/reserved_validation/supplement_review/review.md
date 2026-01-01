# Independent reserved-arm supplement audit

**Arithmetic verdict: PASS for the exact saved reserved evaluation and diagnostics. Scientific outcome: the frozen coverage and overall prediction criteria fail.** These are different conclusions. No calculation or reporting mismatch was found, and no endpoint, exclusion, threshold or model was changed.

The independent verifier reads only the authorized reserved evaluation, reserved diagnostics, frozen sequence-only split, protocol and diagnostic source hash. It does not import the author implementation or open any raw phenotype source. It reconstructs losses from saved predictions and observed targets with exact rational arithmetic and explicit case weights. Graph classifications are checked for complete retention, internal label agreement and correct aggregation; their underlying graph calculations belong to the separate primary review. This is a separately authored review with inherited project context, not a blinded external review.

## Results that constrain the interpretation

| Check | Verified result | Consequence |
|---|---|---|
| Selected inventory | All 24 squares, their 48 designated diagonal identities and all 16 related groups retained. | No missing group or favorable subset selection. |
| Descriptive eligibility | Six squares, 12 diagonals and six groups are eligible. Ten group fractions are null, not zero. | Both requirements, at least ten squares and eight groups, fail. The arm cannot be described as passing validation. |
| Prediction dropout | Seventeen squares lack required known-corner measurements. Seven reach the prediction stage; one is unscorable after unmasking, leaving six scored cases, six distinct targets and six groups. | Missing outcomes are not biological negative cases. All 24 prediction statuses remain in the report. |
| Point restrictions | Two joint-only obstructions out of 12 eligible diagonals. Pooled and equal-group fractions both equal 1/6 because each contributing group has one square. | A selected-cohort descriptive result, not population prevalence. |
| Scenario and range results | One diagonal is joint-only in all 12 dependent deposited-series scenarios. The same diagonal has the fixed-threshold empirical-range certificate. Both fractions equal 1/12. | The arm-wide coverage failure remains. Deposited-series consistency is not independent assay replication, and an empirical range is not a confidence region. |
| Target balancing | Case-balanced and target-balanced summaries are exactly equal as rational values for the full arm and all 16 deletions. | With one scored target per contributing group, there is no distinct weighting experiment. |
| Deletions | All 16 deletions fail coverage and the overall screen. Ten empty-group deletions correctly leave the metrics unchanged. The 127 extrema paths have all 16 values defined; null per-group values remain explicitly retained. | No deletion can rescue the full-data result. These are influence ranges, not confidence intervals. |
| Components | All 12 primary component points are available: nine connected, two joint-only and one single-trait disconnection. Both local obstructions remain disconnected under optimistic observed-component completion. | No observed component bypass, but unmeasured synonymous detours are not excluded. |
| Sensitivities | All 144 named scenario records, 28 separate field/tolerance count tables and six zero flip counts agree. Range tables at each tolerance contain eight connected certificates, one joint-only certificate and three inconclusive cases. | No scenario pooling or artificial replication. No classification/status tolerance flip does not establish calibrated measurement robustness. |

The predicted-but-unscorable case is square 11, group 359. Its masked target `359_21221222_1` has missing `growth.R60one`, `growth.R60two` and `growth.R60three` entries. The saved known-stage missingness reasons for the other 17 squares concern only known vertices. This audit checks the saved records and their consistency; it does not independently re-read the raw archive or the prediction seal.

## Prediction benchmarks and influence

The full-arm group-balanced errors are:

| Frozen method | Reporter-output MAE | Growth MAE | Reporter-output RMSE | Growth RMSE |
|---|---:|---:|---:|---:|
| Additive | 4.9077496492 | 0.01817211092 | 10.0334519411 | 0.02913942505 |
| Adjacent mean | 4.5091593094 | 0.01602234453 | 5.9454529070 | 0.02174054785 |
| Three-corner mean | 7.3026184671 | 0.01987771903 | 10.3235721100 | 0.02642655393 |

The additive predictor worsens both MAEs against adjacent mean, by 8.84% and 13.42%. Against three-corner mean it improves reporter-output MAE by 32.79%, but growth improvement is only 8.58%, below the frozen 10% requirement. Thus three of the four continuous benchmark comparisons fail, in addition to coverage. The additive classifier gets both positive and all four negative cases correct. Both simple baselines miss the two positives, giving group-balanced classification errors of 1/3. These six cases do not justify a broad classification-performance claim.

Omitting group 277 makes all four continuous comparisons pass on the remaining five cases, while coverage and the overall screen still fail. Group 277 supplies additive reporter-output absolute error 24.4390936380 and growth absolute error 0.06808850953. This is an influential failure case that must remain in the full analysis. Omitting group 46 removes the only scenario-consistent and range-certified restriction. Both true classes remain present in every single-group deletion, so no balanced-accuracy deletion value is undefined in this arm.

## The local cases and certification scope

Square 6, group 277: endpoints `277_22321122_1` and `277_31222122_1`. Point joint-only obstruction; not consistent across all 12 scenarios; empirical-range status inconclusive. This square is scored.

Square 20, group 46: endpoints `46_13211213_3` and `46_13321213_1`. Point joint-only obstruction; all 12 scenarios agree; empirical-range status `certified_joint_only_obstruction`. This square is scored. The range certificate uses the frozen fixed low-end thresholds and does not imply universal preservation of endpoint-relative joint-only obstruction after arbitrary interval perturbations. Its uncertainty set has no authenticated statistical coverage.

For protein series One, Two, Three, Four and growth series one, two, three, the separate scenario joint-obstruction counts are `(3,3,1)`, `(3,2,2)`, `(2,3,1)`, `(3,3,2)`, each with denominator 12. Every local scenario obstruction has optimistic observed-component disconnection. These 12 sensitivity fields are dependent and are not 12 biological replicates.

## Verification and adaptation record

The audit checked 7,208 numerical leaves, 4,437 integer leaves and 12,623 categorical leaves. Maximum absolute numerical discrepancy was `1.7763568394002505e-15`, below the unchanged arithmetic comparison tolerance of absolute and relative `1e-12`. These are reproducibility tolerances, not biological thresholds.

The development verifier remains byte-identical at SHA256 `4db97846dba4e165a6cb63861ef369492b74dd531c726f8ef4cc65eee72cfd30`. The reserved copy changes only arm paths/hashes, the frozen arm inventory size from 23 to 24, review labels and additional consistency checks. All 13 non-entry-point functions, including every scientific summary and weighting function, are identical by Python AST comparison. The extra checks compare coverage, dropout and primary evaluation summaries against the same independent arithmetic. `adaptation.diff` and `adaptation.json` record these changes. No scientific criterion deviation was made.

To inspect or replay, use `verify_reserved.py`. It assumes the recorded workspace-relative repository layout, uses the standard library and writes a fresh `independent_results.json`. Existing output is not overwritten. This bounded numerical replay is not a complete raw-source reproduction claim. The full input/output hashes and preservation checks are in `receipt.json`.

Pinned evaluation SHA256: `0f36e4d762b1e80e87859e700f4da2655c6ecf99c29e0560d9d6dcd9546a5521`.

Pinned diagnostics SHA256: `27f0e06eef72012c4d9b2d4d62ce19804359d1678bd26f488f3fd0978d4a60cb`.

Pinned author diagnostic source SHA256: `5e8c1a5638d0d2746f89f12cf26d209f11930f490df99c16a79eda3cf2628850`.

Protocol SHA256: `0f51ccc808559ea69df3bdf294520d81babb67bf14ca92ce2abb6fcc665760ac`; split SHA256: `394a0c1aa23785c290b31c147d199f406cf46b269daa9159a98c419398b05a4e`.

The strongest permitted conclusion is a correctly computed local restriction within this author-preprocessed, incompletely observed library, alongside an inconclusive arm-wide coverage test and a failed frozen prediction screen. The review does not establish RNA mechanism, epistasis, evolutionary constraint, independent-study transfer or a biological null result.
