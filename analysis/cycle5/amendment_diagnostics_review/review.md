# Independent review of amendment diagnostics

Verdict: **BLOCKED for one incomplete amendment requirement**, at diagnostics SHA256 `9e0307ad72b6dd0be771dee9d886bb87e63d09ed85a7ceea5409b2e3f90f897d`. The inspected numerical weighting, deletion summaries and provenance-hold behavior pass the bounded checks below. The missing component/bypass count summaries must be added before claiming that Amendment01 is implemented completely. This is a separately authored review with inherited project context, not blinded review.

The amendment inspected is `cambray_pilot/amendment_01_draft.md`, SHA256 `92c5028dff57df644ad18479f79b4be2e604b7bbe7798e30cdd3dfe5a7fa94d0`. The author tests are SHA256 `07c54d26dcae20a6f2a436e04fb4656452722e4cb19a46bc92346bd341baf5c8`. The staged output-schema source is SHA256 `5fbdf1c0da2dcf9e53b71199b124c576945d432da5b33bbc4240c1b9bfa21ef5`. All four hashes were checked before and after review and remained unchanged.

## Blocking finding

Amendment01 line32 explicitly requires: “Count measured bypasses, certified component results and inconclusive component results separately.” `diagnostics.py` lines228–272 produces per-diagonal records and tolerance `flip_counts`, but its return at line271 contains only `eligible_diagonals`, `flip_counts`, `records` and `scope`. No count of bypass statuses, certified component categories or inconclusive component results is produced. `local_summary`, lines120–146, counts local-square findings only. `analyze`, lines357–363, labels the supplement complete despite this omission.

The retained records make the missing summaries recoverable, but retaining rows is not the same as computing the requested count summaries. This is a reporting-completeness blocker, not a discovered error in the current loss formulas or underlying graph computation. No biological input is needed to correct it.

Add explicit, denominator-bearing aggregate counts alongside the existing records. Keep the following contexts separate because their thresholds and meanings differ:

1. Primary point field at zero tolerance and alternative tolerance: local classifications; available component point classifications plus unavailable point count; component-bound statuses; bypass status conditional on the local joint-only obstruction, with all nonqualifying diagonals counted separately.
2. Each of the12 deposited-series scenarios, separately at both tolerances: the same point/bound/bypass counts. Do not combine these as independent observations. An optional total across scenarios must state that its denominator is dependent diagonal-scenario records.
3. Empirical-range field at both tolerances: local and component certificate-status counts, including inconclusive cases and zero counts. A range-based connection uses a different threshold and premise from the primary point field; do not label it an additional measured primary bypass.

Retain the existing point/scenario records, comparison flags and full-arm primary status. Label current `flip_counts` as classification/status changes, since the implementation compares labels rather than every individual boolean acceptability flag. No outcome-based category selection or successful-case filtering is needed.

## Exact synthetic expectations for the correction

Using the author's existing invented16-group fixture, there are34 eligible diagonals. An independently written counter over the fixture gives:

| Quantity at zero tolerance | Count |
|---|---:|
| Component bounds certified jointly connected |12|
| Component bounds certified joint-only obstructed |10|
| Component bounds inconclusive |12|
| Complete component point classifications jointly connected |12|
| Complete component point classifications joint-only obstructed |10|
| Component point unavailable |12|
| Measured/definite observed-component bypass of a local joint-only obstruction |2|
| Optimistic component disconnection for a local joint-only obstruction |1|
| Component bypass unresolved for a local joint-only obstruction |2|
| Diagonals without a local joint-only obstruction |29|

At alternative point tolerance, bypass counts become2 definite,1 optimistic disconnection,1 unresolved and30 without a local joint-only obstruction. The component bound counts stay12/10/12. Each of the12 scenario fields has these same zero- and alternative-tolerance expectations in this fixture; keep all12 named rows. Both local and component range statuses have2 certified joint-only obstructions and32 certified connections at both tolerances. These are invented fixture values, not measured findings. `expected_component_count_fixture.json` provides the exact named dictionaries for a delta test.

## Checks that passed

The14 author tests were executed from exact byte copies in this review's directory, with temporary inputs redirected to the owned directory. All passed. No actual evaluation or phenotype ZIP was opened. The executable independent reviewer is `verify_independently.py`; it imports only the copied candidate and invented fixture generator, not source outcomes.

The independent check changes the toy signed errors to include cancellation among a duplicated target's predictions. It then constructs exact rational per-case weights directly, rather than calling the candidate's nested averaging helpers. Across the full arm and all16 numeric-order deletions,408 MAE/MSE comparisons agree for both weighting methods, all three predictors and both traits. The hand-calculated additive protein values are case-balanced MAE94/45 and MSE149/30, versus target-balanced MAE43/20 and MSE219/40. The difference verifies that losses are averaged without silently creating an ensemble predictor. Shared targets retain two different endpoint labels.

The equal-contributing-group joint-only fraction is1/18, distinct from the pooled5/34. Empty groups retain null fractions and disclosed denominators. With every synthetic square unavailable, all16 deletion MAEs remain null: minimum and maximum are null, defined rows0, undefined rows16. A missing prediction-record row is refused instead of silently disappearing. Existing one-class and zero-baseline cases fail the predictive screen. The10% hurdle remains exact; numerical1e-12 identity tolerance is not reused as a scientific success tolerance.

The CLI checks the supplied evaluation hash, arm, split provenance hash and held staged-runner schema hash. It refuses existing output, disabled assertions, duplicate JSON keys, nonfinite constants, foreign-arm or missing square inventories and corrupted stored losses. The genuine provenance-hold path was tested with the numerical analyzer replaced by an exception: it was not called, and no scientific summary was produced. Exact evaluation-hash authorization is an investigator-controlled input contract, not a security sandbox against an operator intentionally supplying an unauthorized file.

## Scope and stopping point

This review verifies finite synthetic behavior and source semantics only. It does not recompute connectivity from real measurements or authenticate biological uncertainty, preparation independence or source preprocessing. The mathematical claim is limited to correct arithmetic and preservation of the specified units. It neither proves a biological finding nor tests the actual outcomes.

The review stops at the explicit count-summary blocker. After a narrow correction, review the diff and the expected count fixtures, rerun the held numerical tests, and pin the corrected hash. Existing receipts and this failed-completeness verdict should remain available. No author source, protocol, manuscript, repository or publication was changed by this review.
