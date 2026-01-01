# Supplemental post-unseal diagnostics, synthetic-only handoff

`diagnostics.py` is a portable standard-library script. It has not been run on any real evaluation, phenotype archive or original count matrix. Root must authorize any later real execution after the relevant arm is unsealed. It does not read an archive or any file named by an evaluation payload. Its only scientific inputs are the two explicitly provided JSON files.

```text
python -I -B diagnostics.py --evaluation PATH_TO_AUTHORIZED_EVALUATION --evaluation-sha256 EXACT_SHA256 --split PATH_TO_SPLIT --allow-arm development --output NEW_DIRECTORY
```

Use `reserved_validation` only after root authorizes that stage. Output must not exist. The script rejects `-O` and `-OO` before input access or output creation, a mismatched evaluation hash, the wrong arm, a split hash different from the evaluation provenance, duplicate JSON keys, nonfinite JSON constants, foreign-arm or missing square records, changed group/seed identities, cross-arm shared identities and inconsistent stored losses/labels. The author runner SHA256 is fixed to the reviewed schema revision `5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5`. A successor runner requires an explicit schema review and recorded adaptation, not a silent hash update.

If `interpretation_permitted` is not true, the script writes only a provenance-hold result and input hashes. It does not require or access descriptive/scored subtrees in that path. This permits a genuine held evaluation with no scientific summaries. A supplied exact evaluation hash is root's byte-level input authorization, not a self-issued permission to unseal another arm. These checks are provenance checks, not a security sandbox against an operator deliberately supplying an unauthorized file and its hash.

## Definitions preserved

- Raw descriptive counts use every eligible diagonal. Equal-group fractions average the per-group fractions over groups with eligible diagonals. All 16 selected groups appear explicitly; a zero denominator gives null, not zero. Contributing and selected group counts are separate.
- Primary prediction replay averages case losses within each related group, then groups equally. RMSE is the square root of group-balanced MSE. The target-balanced supplement first averages losses within `(group,target)`, then targets within group, then groups. It does not average predicted phenotype values, vote on labels or merge endpoint-specific truth labels. Repeated targets remain separate classified endpoint comparisons.
- Loss arithmetic is checked against the stored predicted/observed values using fixed absolute and relative tolerances of 1e-12. This only checks numerical identities and is never a biological threshold or a changed primary hurdle. Classification labels are matched to the preserved description for the prediction diagonal. Shared observed targets must match exactly across cases.
- All 16 leave-one-group-out rows are returned in numeric group order, including omission of an originally empty group. Each reports complete retained-case, contributing-group and distinct-target counts, both-class eligibility, continuous/classification losses, additive-minus-baseline MAE differences and the existing coverage/10%/classification hurdles. No model is fitted and no omission is selected. The target-balanced version uses the same case/group eligibility rules; it does not invent a distinct-target coverage threshold.
- A separate all-row minimum/maximum table covers every defined scalar leave-out summary, with defined/undefined row counts. Undefined rows are retained and counted. These are deterministic sensitivity ranges, not confidence intervals, a bootstrap distribution or favorable alternatives to the full-arm result.
- Every eligible descriptive diagonal retains the point and range statuses, observed-component bounds, native shortfalls, all12 preparation scenarios, numerical tolerances and alternatives. A definitely acceptable observed-component connection is a bypass witness status; an optimistic disconnection remains limited to that observed graph. Unknown bypasses are labeled unresolved. No conclusion is extended to unmeasured synonymous genotypes.
- Flip counts compare zero versus alternative classification/status labels: local point classification, the component `(point classification,bound status)` pair, local/component range status, and the corresponding 12-scenario labels. Full flags remain in the record even when a status label does not change. These counts do not count the threshold surface as independent observations.

The original evaluation remains unchanged. Stored primary summaries are not overwritten or relabeled. This supplement reconstructs new summaries from preserved row-level results; it does not independently re-evaluate graph connectivity from raw measurements. The author's existing primary output and its implementation review remain the source for those graph classifications.

## Synthetic evidence

```text
python -I -B test_diagnostics.py
```

All 14 tests passed on CPython 3.12.14. The fixtures have 16 groups, 18 selected squares, 17 scored cases, one zero-eligible group, overlapping targets with different endpoint labels, missing-target/known-input failure statuses and 12 scenarios per eligible diagonal. No real IDs or phenotype values are used. Hand calculations give primary group-balanced additive protein MAE `97/45`, target-balanced MAE `67/30`, primary RMSE `sqrt(259/45)` and target-balanced RMSE `sqrt(6.6)`. The descriptive raw joint-only fraction is `5/34`, while the equal-contributing-group fraction is `1/18`. They intentionally differ.

Tests also cover one-class discrimination refusal, a zero-error comparator, nine-case coverage failure, no scored cases, all numeric-order leave-outs and their ranges, field corruption, missing scenarios, cross-arm injection/refusal, hash/schema/overwrite refusal, input immutability, provenance hold and optimization refusal without an output. An initial 12-test pass was followed by two additional explicit one-class/coverage tests and defensive flag/type checks, then a 14-test pass. There were no failed test runs or biological executions.

The source schema, split/mask and owned code/tests are hash-recorded in `handoff_manifest.json`; `synthetic_receipt.json` records test scope. The code is held for root review. No protocol, staged reader, manuscript, outcome, source report or repository file was edited by this task.
