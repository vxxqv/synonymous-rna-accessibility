# Independent Cambray preflight review

Status: PASS at the held source hashes, conditional on the primary agent's explicit protocol/source freeze and staged execution authorization. No blocking implementation defect was found. This document does not itself authorize observed phenotype access.

## Scope and independence

The independent evaluator uses direct Boolean transitive closure, independently constructed weak threshold cells, explicit interval bounds and an independent equal-group metric implementation. `reference.py` does not import the primary implementation and performs no file reads. A separate subprocess probe imports the held primary API only after independent expected results are calculated. No phenotype archive, deposited phenotype table, observed arm file or observed result has been opened in this review.

The exact reviewed protocol is SHA256 `9becbd723fac30ac12ad9e95d5fb80cb69b5bd93c99210668e0479a24879adde`. The initial independent freeze is retained. The v2 freeze records the pre-outcome update to the protocol's explicit floating comparator `f >= tau - delta`, rather than the real-equivalent but floating-distinct `f + delta >= tau`.

## Completed checks

Twelve independent unit tests pass, including additive antagonism, weak boundaries, all four withheld-corner orientations, strict three-known-corner interfaces, nonmonotone bypasses, incomplete fields, nonnested preparation and interval criteria, group weighting, zero-error baselines and single-class classification failure.

The held primary `graph_rules.py`, SHA256 `7f901518cfc196de5d6d24c10f6f1b336255063ca8bd22d544accb6ef7b940f3`, agrees on 14,018 synthetic API invocations: 13,122 point fields, 625 interval fields, 81 additive projections, 108 local predictions, 81 critical-cell geometries and one complete twelve-preparation bundle. Geometry checks include 2,916 threshold-boundary probes. Input hashes match before and after the comparison. This is finite implementation evidence, not a proof of arbitrary-input correctness.

An independently handwritten standard genetic code and direct Hamming graph reconstruction agree with the sequence-only inventory: 289 vertices, 43 observed components, 299 edges and 47 squares. All 244,000 sequence metadata records were checked under an exact sequence-only header whitelist. Independent peptide comparisons and SHA partition/mask reconstruction agree. Development has 23 descriptive squares, 22 union-mask IDs and 22 eligible prediction squares with 21 distinct targets in 16 groups. Validation has 24 squares and 24 distinct masked targets in 16 groups. Related-group dependence and the repeated development target are retained.

Full machine-readable receipts are `graph_comparison_v1/comparison.json`, `metadata_check_v1/metadata_check.json` and `reference_freeze_v2.json`.

## Staged-reader review

The held runner is `cambray_pilot_preflight/staged_analysis.py`, SHA256 `5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5`. Its author test source is SHA256 `8e7684aeb42fbfb8a710f66e9c1fefba569b9d3eafa0595f7f730e55ffe6bf37`. All 13 author tests pass. An independently constructed twelve-square, eight-group fixture also passes, with all four withheld-corner orientations and a partially measured component-only vertex. The independent receipt is `staged_check_v1/receipt.json`.

The following code transfers were inspected directly:

| Location in held runner | Review conclusion |
|---|---|
| Lines27-28 | Explicit refusal under `-O` and `-OO` occurs before imports, input access and output creation. Both flags were also tested in subprocesses. |
| Lines65-90 | Every declared input hash and executing code path is checked. Reserved validation additionally requires a pinned development report, protocol and unchanged-code decision. This is an external authorization gate, not independent authentication of historical exposure. |
| Lines93-150 | Row identity is read first; excluded rows' measurement bytes are skipped. Only explicitly selected columns of allowed IDs are decoded. Duplicate source IDs, absent selected IDs, ambiguous CSV members and source record-count disagreement fail. The format assumes one physical row per record, consistent with the authenticated pinned source; this is not a general CSV parser. |
| Lines153-167,193-225 | The prediction ID set contains only the three visible corners of fixed exactly-one-union-mask squares. Author means, component maxima and withheld finiteness are absent from that reader request. Predictions and known inputs are sealed with all base input hashes. Metadata-ineligible squares retain explicit prediction exclusions. |
| Lines207-216,309-313 | The classification endpoints are the two corners adjacent to the withheld intermediate. All four possible target orientations agree with the independent evaluator. Negative additive predictions are retained. |
| Lines229-239,350-356 | Both numerical-payload hashes and base input hashes are rechecked before evaluation calls the outcome reader. Tampering with either sealed numerical file was rejected while the reader was replaced by a poison function. |
| Lines170-187,247-287 | A trait mean or range is unknown if one of that trait's required preparations is missing; the other trait remains available. Point, preparation and range tolerances are derived separately over the component. Known failures remain failures in optimistic joint bounds. |
| Lines297-347 | All cases and failure reasons are retained. Metrics average within related groups, then equally across groups. Both traits and both simple baselines enter the ten-percent hurdle. Zero baseline MAE cannot pass; one-class data cannot pass discrimination. |
| Lines359-380 | Author-mean checks occur only after sealing and use exactly the declared maximum absolute/relative tolerance. Missing author means or larger disagreement hold interpretation, without silently changing means or thresholds. |
| Lines389-433 | Both diagonals and all twelve preparation scenarios are retained. Whole-component bypass and unknown-state bounds are reported. The per-arm milestone is not the final two-arm decision; the primary agent must combine both separately sealed arm reports under the protocol. |

The independent fixture matches all three methods' predictions, both-trait group-balanced MAE and RMSE, classification errors, confusion tables, balanced accuracy and the complete screen. It also matches 24 diagonal point/range decisions and 288 whole-preparation flags. Component-derived point tolerances were checked explicitly. Invalid UTF-8 numeric bytes in all masked targets, author means, a component-only vertex and an unrelated row do not prevent prediction, demonstrating that those bytes are not text-decoded in that phase. The normal fresh-output-only refusal and both optimization refusals pass. These tests use invented data only.

Missing-target behavior, missing-known-input behavior, author-mean holds, absent development decisions, source hash refusal, quoted-comma parsing, nonfinite values and zero-error/one-class failures also pass the retained author synthetic suite. The independent graph oracle separately checks the mathematical quantities used by that suite.

## Qualifications and release conditions

The historical standalone masking and metadata drivers contain assertion-based guards without their own explicit optimization refusal. Their normal-run retained results were independently checked. The reviewed staged entrypoint resolves this concern for its own execution by rejecting optimization before importing those helpers, then checking their hashes and the pinned mask. Do not treat an unguarded optimized standalone metadata-driver replay as certified by this review.

Metadata is trusted only after its independent sequence-only verification and exact root freeze. The runner is not a security sandbox and cannot prove that no human or separate process previously saw a target. Hash checks also do not eliminate file-system races with a concurrent writer; inputs must be held unchanged during execution. Reserved execution is contingent on the root's separately pinned development report and unchanged-procedure decision. Validation must not be opened merely because this implementation review passed.

No biological calibration, independent-preparation inference or population generalization follows from these synthetic and metadata checks. In particular, twelve preparation combinations are dependent sensitivity scenarios, and neither threshold geometry nor additive antagonism alone establishes a biological mechanism. The finite tests do not replace the exact graph arguments or the protocol's uncertainty qualifications. The final verdict is a bounded implementation preflight pass, with zero observed phenotype sources opened.

## Reproduction

Run from the workspace root with the existing interpreter and bytecode writing disabled:

```text
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_independent_review/test_reference.py
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_independent_review/compare_synthetic.py --primary work/joint_function_cycle5/cambray_pilot_preflight/graph_rules.py --out NEW_GRAPH_OUTPUT
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_independent_review/check_metadata.py --preflight work/joint_function_cycle5/cambray_pilot_preflight --eligibility work/joint_function_cycle5/cambray_eligibility --out NEW_METADATA_OUTPUT
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_pilot_preflight/synthetic_reader_tests.py
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_independent_review/check_staged.py --primary work/joint_function_cycle5/cambray_pilot_preflight/staged_analysis.py --out NEW_STAGED_OUTPUT
```

Every output argument must be new. The independent tests have no observed-phenotype input option. `reference_freeze.json`, `reference_freeze_v2.json` and all earlier receipts remain unchanged as history.
