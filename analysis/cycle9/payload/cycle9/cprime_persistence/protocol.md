# Exploratory persistence and concentration of marker sensitivity

## Question and prior

Among the complete already exposed 30 intended labels, which assignment-robust marker-dispersion conclusions persist in both cell contexts, and are they confined to one declared exon? The original cPRIME paper already demonstrates selected marker effects and cell-context differences. Cycle 8 already reports all 60 individual profiles and eight assignment-robust profiles. This is a complete cross-tabulation of those accepted results, not a new independent biological study or discovery of marker non-neutrality.

The potentially useful increment is a precise distinction between robustness within one context and transport across two measured contexts, including denominator and exon concentration. A concentrated or nonpersistent result would narrow the proposed central argument. No inferred marker identity, phenotype-fitted feature, or predictive model is used.

## Fixed inputs and eligible units

Use the complete Cycle 8 cprime_followup/results_v1/result.json, SHA256 851bbf75b66dd0186d0b99d7536409773de1afdae5c3ea83c57563655d55308e, plus its original B1 source SHA eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf and metadata-only intended-label/exon mapping. Authenticate all hashes and joins before reporting. Retain all 30 intended labels, both K562 and KBM7, all five declared exons and all 60 profiles. Any ambiguous intended-to-exon join is retained as unknown and stops exon interpretation; no phenotype-based exclusion.

Outcomes and individual-profile findings were previously examined. This protocol is written before the new aggregate cross-tabulation, but is explicitly exploratory and cannot create prospective validation. Dependence includes shared markers, overlapping local sequences, targets within exons and two cell contexts from one source. No confidence interval, prevalence estimate, p-value, bootstrap replication or train/test split is justified here.

## Estimands and comparisons

For each target and cell retain exact A=max L/min U, B=max distinct U_i/L_j, assignment count, k_min and k_max from the verified finite assignment calculation. Use the unchanged twofold display criterion. Classify:

- universal-three: k_min=3;
- universal-any: k_min>=1;
- assignment-dependent count: k_min<k_max;
- no-twofold: k_max=0;
- continuous equality exclusion: A>1;
- continuous guaranteed-twofold: A>=2.

These flags overlap. For every flag report the full paired-cell 2x2 table, all target-level records and all five exon strata with their denominators. The intersection across both cells is a cross-context persistence statement for the finite declared cohort. It is not equivalent to shared direction, marker ordering, additive representability or mechanism. The union and exclusive-cell cases are retained.

Report each exon’s contribution to the universal-three profiles and to the both-cell target intersection. Compute all five leave-one-exon-out descriptive counts without refitting or choosing the favorable omission. These are concentration diagnostics, not uncertainty estimates or independent replication.

## Verification and decision rules

Before actual summary access, verify synthetic cross-tabs covering all four paired-cell states, exact equality at two, missing/duplicate cells and inconsistent exon maps. Independent reviewer reconstructs all flags from exact A and explicit assignment folds where available rather than accepting stored summary booleans alone. Main-text admission requires complete denominator reconciliation and no discrepancy. If new grouping is ambiguous or the interpretation merely restates the existing counts, retain the cross-table as a supplementary description and do not claim a new biological mechanism.

No expansion into motif association, location trend or sequence model is authorized by this protocol. Such analyses would need authenticated features, a separate question and a credible design capable of distinguishing alternatives.
