# Metadata overlap and dependence audit

The fixed split passes the recorded design-overlap check. All 552 cross-arm square pairs have zero shared square vertices, exact nucleotide-mutation edges, original masked targets, global masked IDs within squares, peptides, design seeds, related seed groups, observed components or component vertices. This supports a separation of the analyzed design sets. It does not create independent assay batches, libraries or biological experiments.

No phenotype ZIP was opened. The audit reads only five hash-pinned metadata inputs and checks all 1,081 pairs of the original 47 squares. It retains all 94 descriptive diagonal pairs. It independently reconstructs the four mutation edges from each square's sequences, checks those edges against the full saved edge inventory, verifies each selected observed component is complete relative to the saved edges, and rechecks identical translated peptides with Biopython. It does not repeat an all-library relatedness search. That earlier analysis-set relatedness result is an existing frozen input, independently checked in the separate preflight review.

| Recorded design quantity | Development | Reserved |
|---|---:|---:|
| Squares / diagonal pairs |23 /46|24 /48|
| Related seed groups |16|16|
| Distinct square vertices |88|93|
| All observed-component vertices |133|156|
| Distinct square mutation edges |90|95|
| Prediction cases before outcome missingness |22|24|
| Distinct prediction targets |21|24|
| Pairs sharing at least one square vertex |2|2|
| Pairs sharing an exact mutation edge |2|1|
| Pairs sharing a peptide |2|3|
| Pairs in the same seed/group |8|8|
| Pairs in the same observed component |2|2|

All pair counts above are unordered square pairs, not biological replicate counts. `square_pair_relations.csv` includes every pair, including all negative relations. `square_degrees.csv` gives each square's degree for every relation. `group_inventory.csv` and `prediction_target_multiplicity.csv` expose all grouping and target counts.

The four vertex-overlap pairs are development squares36/37 and44/45, and reserved squares30/31 and34/35. Development36/37 share the target `227_23132123_1`; both are in group227. Development44/45 share two vertices, one of which is also masked in another square; the already frozen mask excludes square45 from prediction only. Reserved34/35 share one vertex but no mutation edge. Reserved15/16 have the same peptide without shared vertices or an observed component. These patterns were not used to change the split or choose outcomes.

## Concentration and dependence limits

Development has ten groups with one descriptive square, five with two and one with three. The largest group's pooled-square share is3/23, or13.0435%. Reserved has eight groups with one square and eight with two; its largest pooled-square share is2/24, or8.3333%. The respective sums of squared group case shares are39/529 and40/576. These are structural concentration summaries, not effective sample sizes.

For the metadata prediction cases, development has ten groups with one case and six with two; reserved has eight with one and eight with two. The current equal-group prediction rule gives every represented group weight1/16, or6.25%, before outcome missingness. Equal weighting does not make groups independent. The primary reporter-output and growth preparations, shared preprocessing and common study are the same across both arms.

Source-level dependencies are retained from `cambray_eligibility/report.md` and `provenance.json`: four protein preparations with shared normalization; three growth competition lines sharing an initial culture and a reference; no authenticated matching of protein preparations to growth preparations; a common reporter backbone and experimental library. The12 preparation combinations are robustness scenarios and are not12 independent experiments. The numerical covariance induced by preparation, normalization and shared reference is not recoverable from sequence metadata and is not estimated here. No source-independent assay validation is claimed.

## Fixed sensitivity choices proposed before outcomes

These are reporting proposals for root's explicit protocol amendment. They do not modify the held runner, current success criterion or split. Missing cases remain in the coverage table with their reasons. An unavailable loss is never assigned zero. All denominators are recomputed from the disclosed scorable cases and accompanied by the original metadata denominator.

1. **Equal related-group loss, retained as primary.** For each trait and each method, average each scorable case's absolute error within its related seed group, then average the represented groups equally. For RMSE, use squared errors at both averaging steps and take one square root afterward. Apply the same hierarchy to binary misclassification error. Keep pooled confusion counts as a separate dependent-case description. For the descriptive endpoint, retain the protocol's raw numerator/denominator, and additionally report each group's diagonal fraction and their equal-group mean. A group with no eligible cases contributes an explicit unavailable row, not a favorable zero.

2. **Target-balanced sensitivity without dropping cases.** Within each group, first average the losses of every case sharing a target ID. Average those target-specific losses equally, then average groups equally. Do not average predicted values before taking losses, which would define a different ensemble predictor. For classification, retain each case's own endpoint pair and truth label; a shared target does not imply a unique binary classification. Report this as target-balanced case loss, plus every original case and its weight. The design has only one duplicated target, and its group contains exactly those two cases. Every other group has one target per case. Consequently target-balanced loss is algebraically identical to current equal-group case loss for this fixed metadata, including any subset left by missingness. It will not be a second independent confirmation. The distinct-target denominator remains21 rather than22 in development.

3. **Leave-one-related-group-out sensitivity.** Enumerate each of the16 group deletions in fixed ascending numeric group order, in each arm. For every method and trait, recalculate the same equal-group and target-balanced errors, all comparison differences and complete descriptive numerator/denominator. Keep all16 rows and disclose the minimum and maximum over them together. Do not select a favorable deletion, change a target, refit, retune or substitute a deletion result for the primary screen. Report retained squares, targets, groups, both-class coverage and any operational coverage failure for every deletion. Such deletion sensitivity assesses concentration within this design set, not a confidence interval or future-study performance. The full-data criterion remains the criterion unless root explicitly freezes an additional pre-outcome rule.

The deduplication identity follows directly: a group either has distinct targets, so averaging losses by target equals averaging cases, or consists entirely of cases for one target, so the target average also equals the case average. The fixed design has no group mixing duplicated and singleton targets. Removing cases cannot create such a mixture. This is a mathematical accounting fact; no measured error is needed.

## Reproduction and limits

Use the existing CPython3.12.14 and Biopython1.88 environment. The script rejects disabled assertions, verifies all five input hashes, has no phenotype reader and refuses an existing output directory. Its input paths are explicit:

```text
python audit_overlap.py --split SPLIT_MANIFEST.json --mask PREDICTION_MASK_V2.json --components COMPONENTS.jsonl --sequences SEQUENCE_METADATA.tsv --edges SYNONYMOUS_EDGES.tsv --output NEW_DIRECTORY
```

`metadata_v1/receipt.json` records the code/input/output hashes. A newly discovered cross-arm relation would fail this audit before writing an accepted result and require pre-outcome handling. This audit found none. The already known shared assay dependence remains a limitation, not a reason to mislabel the split as an independent experiment. Outcome coverage, class balance, empirical margins and actual predictive performance remain unexamined.
