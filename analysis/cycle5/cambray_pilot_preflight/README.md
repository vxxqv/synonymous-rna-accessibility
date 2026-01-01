# Cambray metadata and synthetic preflight

Status: ready for protocol freeze and separate-person review. No new phenotype value has been decoded, modeled or classified. This package reads only the existing sequence-only inventories. It has no phenotype archive reader, fitting routine or real-outcome execution command.

## Fixed split

The root-approved clarification compares only unique peptides in the47 selected squares and their observed synonymous components. Every seed1 peptide is used for the prior-exposure relatedness exclusion. The denominator is32 amino acids; at least70% identity means at least23 matching positions. This is analysis-set relatedness, not universal phylogenetic or assay independence.

All851 cross-seed selected-peptide comparisons and117,726 comparisons against seed1 peptides have at most7 matching positions. No selected seeds join and none is excluded. The original individual-seed hash rule therefore applies, with the first16 seeds in SHA256 order assigned to development and the remaining16 reserved. No coverage requirement or success threshold was changed.

| Arm | Seeds | Squares | Diagonal pairs | Peptides | Observed components | All component vertices |
|---|---:|---:|---:|---:|---:|---:|
| Development |16|23|46|21|21|133|
| Reserved validation |16|24|48|21|22|156|
| Total |32|47|94|42|43|289|

Every original square and both diagonals are retained in `split_manifest.json`. The manifest also fixes the withheld corner by the approved lexical UTF-8 hash convention, coordinate order, classification endpoint pair and isolated-component stratum. There is no crossing of shared vertices, peptides, components or design seeds between arms. Four isolated observed squares lie in development and six in the reserved arm. Actual finite phenotype coverage remains unexamined; these metadata counts are only upper bounds on eligibility.

The alternate checker reconstructs components from the edge list, checks every square's actual nucleotide distances and coordinate labels, recomputes hashes independently and uses a vectorized amino-acid representation to corroborate relatedness. It agrees exactly. This is distinct-algorithm verification by the same analyst, not blinded or separate-person review.

## Reusable functions

`graph_rules.py` is input-agnostic. It provides:

- Simple undirected graph validation and reachability.
- Endpoint-relative or explicit-threshold point classifications, with weak inequalities and ineligible endpoint status.
- Pessimistic/optimistic interval certificates and incomplete-vertex bounds. Missing traits are represented by `None`. A known failing trait still excludes a vertex from joint optimistic acceptability.
- Component-relative numerical tolerance values as an explicit sensitivity, with exact-zero tolerance the default.
- Exact critical-value atoms and open threshold cells for two traits. Cell representatives use rational arithmetic; cells carry no probability or statistical weight.
- The four-corner Hadamard additive projection and interaction contrast.
- Three missing-corner predictions accepting exactly three known corners. Passing the held-out target as an additional dictionary entry is rejected. Negative additive predictions remain negative.
- Native-unit shortfall pairs and the fixed corner hash.

Callers must supply complete intended graphs, the documented corner order00/01/10/11, finite known measurements, explicit unknown entries, correct assay conditions and the applicable frozen thresholds. Arbitrary sequences are never treated as measured values. The interval rules are sufficient certificates under the supplied ranges, not probabilistic uncertainty estimates or necessary characterizations under correlated uncertainty. Fixed external thresholds that do not retain the endpoints cannot produce a barrier classification.

## Synthetic verification

The11 tests pass. Their independent reference calculations include all16 intermediate patterns,13,122 exact comparisons for all ternary two-trait square fields and both diagonals,256 explicit realizations of a certified interval example, and500 exact rational additive missing-corner recoveries. Additional cases cover aligned traits, a single-trait disconnection, ties, numerical tolerance, ineligible endpoints, a Q3 detour, unknown detour vertices, known failures mixed with unknown traits, nonfinite and invalid inputs, threshold-cell boundaries, monotone transformations and un-clipped negative predictions.

An independent linear-algebra calculation verifies rank3 of the square additive design and residual dimension1. The three-state design has rank3 and residual dimension0, so its perfect in-sample additive fit cannot identify interactions. The additive antagonism control has zero interaction contrast but joint disconnection. These results verify implementation and identifiability; they are not biological evidence.

Commands from the workspace root, using the existing environment:

```text
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_pilot_preflight/metadata_split.py
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_pilot_preflight/verify_split.py
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_pilot_preflight/test_synthetic.py
```

These commands reproduce preflight metadata and synthetic receipts only. They overwrite their own preflight receipts and must not be substituted for a future guarded biological execution stage. The protocol draft hash in the split manifest records the draft at calculation time; the root must pin the approved final protocol separately.

## Clarifications and remaining limits

The root resolved relatedness scope, group hashing, corner-string serialization and assignment direction before this split. The retained fallback group convention is recorded but unneeded because no selected seeds join.

Two scoring/provenance details were returned for clarification and are deliberately not implemented here: the numerical agreement rule for archived rounded author means, and whether classification-error noninferiority is seed-balanced against both simple prediction baselines. The root owns those final rules. This preflight does not parse measurement columns, evaluate coverage, generate biological predictions, compute biological metrics or authorize phenotype access. Separate-person review and the root's final freeze remain required.
