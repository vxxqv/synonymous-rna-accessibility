# Independent cPRIME persistence review

**Numerical verdict: PASS, exact reconstruction. Release guard: one hardening change requested.** No scientific discrepancy was found. The current author `analyze.py` uses assertions for its source hashes, joins and scientific checks without explicitly refusing optimized Python. `-O` or `-OO` removes those checks. Add an explicit `if not __debug__: raise RuntimeError(...)` before imports or any input/output for reliable replay. This does not require changing results, models, thresholds or the exposed-cohort exploratory status. This review pins the unmodified source; a later guard-only delta needs a separately recorded hash.

The independently authored `verify.py` does not import production functions. It first reconstructs each exact response from its paired B1 RNA and DNA fractions, groups by intended label/cell/state, orders rows by source-row identity, and enumerates relative permutations while fixing only the first marker's common label order. Independent scientific summaries are complete before reading B2 summary fields or the new author result. It then compares all exact folds, marker/source-row joins and output fields.

Coverage: all 396 paired rows, all 66 local states in two cells, all 60 profiles, all 720 assignments and 2,160 recorded folds agree exactly. All 60 exported profile records, six complete paired-cell tables, five exon strata and five leave-one-exon-out summaries agree. No profile or marker was excluded. All pinned inputs were unchanged after verification. The reviewer uses exact rational weak comparisons at two and a strict comparison above one; no floating tolerance is applied.

| Classification | Neither cell | K562 only | KBM7 only | Both cells |
| --- | ---: | ---: | ---: | ---: |
| Universal three recorded twofold dispersions | 24 | 4 | 0 | 2 |
| Universal at least one | 6 | 5 | 3 | 16 |
| Assignment-dependent count | 12 | 5 | 8 | 5 |
| No assignment with any twofold dispersion | 21 | 3 | 2 | 4 |
| Continuous equality exclusion | 16 | 6 | 3 | 5 |
| Continuous guaranteed twofold | 25 | 3 | 0 | 2 |

The denominators are 30 intended labels, not 60 independent experiments. The universal-three intersection comprises `chr17-31225232-A-T` in exon 17 and `chr17-31233126-A-G` in exon 27. The eight qualifying profiles occur as five in exon 17, two in exon 27 and one in exon 39. Exons 25 and 42 contribute none. Removing exon 17 or 27 leaves one both-cell target; removing any other exon leaves two. Every omission leaves at least three universal-three profiles. These are concentration descriptions, not evidence of exon independence or a population enrichment test.

Eighteen independently written synthetic checks cover all four paired-cell quadrants, exact threshold ties and values either side by 1e-20, missing/duplicate cells, inconsistent exon maps, invalid/missing profile rows, positive-domain guards, permutation invariance and preservation of numerator/denominator pairing. Full permutations of all markers agree with the fixed-first-marker quotient with exactly sixfold multiplicity. An exact example with recorded triples (1,4,16) and (2,8,32) has universal-three discrete dispersion even though continuous intervals overlap. This checks that the two uncertainty statements are not silently identified. The independent verifier also refuses `-O` and `-OO` before input access.

The Cartesian interpretation is appropriate for the declared possibility set, not a probability distribution. Anonymous row labels are independently reassigned across markers; authentic preparation keys are not recovered. The both-cell conjunction states that each cell separately satisfies the same dispersion predicate under every permitted within-cell assignment. It does not require, infer or recover cross-cell replicate pairing. If an authenticated assignment set were smaller, universal statements over the current larger set would remain conservative; unsupported additional observations or error allowances need their own analysis.

Persistence of fold dispersion is not preservation of marker direction or rank. For example, two markers with constant recorded triples (1,1,1) and (2,2,2) in one cell, reversed in the other, satisfy universal-three in both cells while reversing marker order. Nor does either result identify an isolated intended effect: marker main effects alone suffice, and unrestricted edit-marker interactions preserve the known isolated-effect shift ambiguity.

This is descriptive exploration of an exposed, selected cohort within one source. It is not independent biological replication, calibrated simultaneous confidence, a prevalence estimate, common-order preservation, additive representability, or a mechanism. The original source already establishes marker non-neutrality. The incremental result is a complete, denominator-preserving concentration and persistence summary of the previously verified finite assignment results. Both-cell cases in two exons narrow the statement that the finding is purely one-cell or one-exon, without supporting broad context transport.

No source workbook was opened, no threshold or normalization changed, no phenotype model was fitted and no shared file was edited. `run_v1/verification.json` records exact source/code pins and all reconstructed summaries; `run_v1/independent_profiles.csv` preserves every profile.
