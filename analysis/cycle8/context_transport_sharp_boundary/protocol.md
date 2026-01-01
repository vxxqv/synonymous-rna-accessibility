# Exploratory D2: exact shared-descriptor transport boundary

Status: proposed protocol, held before any new boundary fit or outcome-derived breakpoint computation. Root review and explicit execution authorization are required. Existing B1/D1 outcomes are exposed; this is exploratory reuse, not new validation. Ownership is this directory only. No workbook, source counts, new genotypes, endpoint, descriptor, response class, tolerance or biological model is introduced.

## Fixed question and inputs

Does existing M0, the shared descriptor-additive latent score with a separate arbitrary nondecreasing response in each cell/exon, become feasible at exactly the previously verified arbitrary-common-order widening boundary? If so, the added descriptor restriction makes no stricter threshold on the specified source design. Failure at that point requires a finite exact search for its own threshold, not a change of model.

Retain all 66 local states, five exons 17,25,27,39,42 and both contexts K562 and KBM7. Retain the full existing M0 matrix, including its original column ordering and exon blocks. M1 is the original context-specific descriptor-additive model. State names denote recorded local allele descriptors, not authenticated full haplotypes or invariant protein products.

Required immutable inputs:

- B1 result `work/research_cycle8/cprime_results/primary_v1/result.json`, SHA256 `eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf`.
- Matrices `work/research_cycle8/context_design_audit/inventory_v1/nested_matrices.json`, SHA256 `665c4d0586310d09868b07f01163875e82708c3441201ab7757d7c897b3ab874`.
- Metadata `work/research_cycle8/context_design_audit/exact_v1/mapping.json`, SHA256 `a2f428f5308cb78a827bfbd23f7d7756b98980e3585f0ccf6a733859c2388436`.
- Exact arbitrary-order boundary result `work/research_cycle8/cprime_followup/results_v1/result.json`, SHA256 `851bbf75b66dd0186d0b99d7536409773de1afdae5c3ea83c57563655d55308e`.
- Existing D1 solver `work/research_cycle8/context_design_audit/nested_model/runner.py`, SHA256 `2a58c52ee80ea234bc1466a207ae5e4f9b020f142088721b239abcf1a5b86cc9`; its frozen protocol, addendum, freeze and actual result are pinned in the companion input manifest.

Decode B1 only through the existing selective parser: state/cell positive interval endpoints, completeness and source-row identity; source-row exon metadata; prior grid statuses. Do not decode RPM, source-score, ratio, median or group-result fields. Inspect the boundary result only for exon boundary/graph/witness records required for exact lower bounds. No source workbook is opened.

## Parameterization and decisive proof

Let `t=gamma^2 >= 1`. Replace symmetric intervals `[L/gamma,gamma U]` by positively rescaled intervals `[L,tU]` in each context. Because response functions are unrestricted nondecreasing and context-specific, this uniform positive output rescaling does not change model feasibility. Unknown intervals remain `(0,infinity)`. Finite endpoints and decisions use exact fractions.

Within a context, `i` precedes `j` only when `t U_i < L_j`. For M0 require `(X_j-X_i) beta > 0`, equivalently the finite unit-margin system after scaling beta. Response functions remain nondecreasing, not strictly monotone. At equality `t U_i=L_j`, that pair is not a strict constraint.

M0 implies an arbitrary common latent order. Thus each previously certified arbitrary-order boundary is a lower bound for the M0 boundary. If an exactly checked M0 witness exists at that bound, it establishes equality of the two thresholds. This is the preferred terminating test, not a grid interpolation.

## Algorithm and bounded fallback

1. Verify every input hash and runtime before parsing; verify all 132 interval records, 396 selected row identities, 66 local-state identities, matrix row order and exon blocks using the existing held contract reader. No source re-estimation.
2. Recover the five existing sharp t values from the pinned arbitrary-order result. Validate every retained lower-bound cycle and boundary order witness against B1 intervals directly. Reuse previously reviewed lower-bound mathematics only after exact edge checks.
3. For each exon, test M0 at its existing boundary using `[L,tU]`. Candidate search may reuse the byte-unchanged existing solver, after its hash is checked. Every result requires an independently reconstructed exact certificate: all rational beta entries, forced pairs, positive gaps and contextwise monotone response completions, or all positive dual weights and exact cancellation. A floating optimizer status is never a result.
4. Verify from the pinned matrix that exon 27's single-context M0 design has row rank equal to its number of states, using rational elimination. Retain pivot rows/columns. Full row rank means every assigned common scalar vector is representable as X beta; an arbitrary-order witness can therefore be lifted by an exact linear solve without scientific fitting. This proves equality for exon 27 under nonempty intervals and this declared design. Do not infer this property merely from a rank label in a report.
5. The other three exons have arbitrary-order boundary one and existing M0 witnesses at one. Validate those witnesses exactly rather than refit unnecessarily. Verify M1 feasibility at one in all five exons. Feasibility remains true at all greater t because intervals are nested. M1's non-rejection is not biological support for additivity.
6. If and only if exon 17 M0 is exactly infeasible at its arbitrary-order boundary, form the finite set of within-context rational breakpoints `L_j/U_i` in `[t_arbitrary,25/16]`, plus the lower bound and `25/16`. Unknown endpoints contribute no finite breakpoint. Use all such points, with exact duplicate removal. A monotone binary search may locate the first feasible point; retain every evaluated point and certificate. No point above the existing feasible `25/16` is permitted, and no new scientific tolerance is used.
7. Certify minimality in that fallback by an exact infeasibility certificate for the immediately preceding open chamber: evaluate its rational midpoint and save a positive dependence. Every edge in that certificate must persist throughout the open chamber and until its right boundary. Check this by its exact edge-loss thresholds. At the first feasible breakpoint, save an exact primal and response completion. Thus the dual excludes all lower t by nesting and the primal proves attainment. If the lower bound itself is feasible no chamber search is needed.

The matrices do not depend on t; the forced comparison set is constant on every half-open chamber between consecutive endpoint-ratio breakpoints. At the right endpoint all exactly tied comparisons are removed. Therefore this finite procedure, with exact certificates on both sides, proves a sharp threshold. It is a standard finite ordinal linear-feasibility algorithm, not a novel generic method.

## Implementation and audit requirements

Runtime is the already used CPython 3.12.14, NumPy 2.5.4, SciPy 1.18.1. The existing rational reconstruction settings may be reused unchanged: denominator cap 1000000 for a candidate primal; floating dual support cutoff 1e-9 only to propose a rational dual, followed by exact elimination, positivity, unit-sum and complete cancellation. Failed reconstruction is unresolved and stops execution. No tolerance is applied to forced orders, ties, endpoint ratios or final decisions.

Before exposed-input execution, seal implementation hashes and run invented tests covering: a rank-full lift; a rank-deficient design with M0 strictly more restrictive than arbitrary order; an edge disappearing exactly at its breakpoint; a feasible no-forced-order case; unknown intervals; dual/primal tampering; row/column permutations; no overwrite; hash drift; and `-O/-OO` refusal. Preserve failed tests and attempts without rewriting earlier evidence.

Every run uses explicit input paths and a fresh output directory, records protocol/code/input/runtime hashes and validates inputs again after computation. No network, Git or source changes occur. On any matrix/identity/certificate/rank contradiction, stop and preserve the result as unresolved rather than revise the model or choose a favorable case. Root receives the complete report and exact certificates for fresh independent review before manuscript admission.

## Required complete output and interpretation

Retain all five exons, their existing arbitrary-order boundaries, M0 boundary decisions, any exact rank lifts, all fallback breakpoints and evaluations, all exact witnesses/duals and response checks, and all M1-at-one checks. Report exact t before bracketed gamma displays. Square-root brackets, if shown, must be verified by rational squaring and are display-only.

The possible conclusion is either exact equality of the declared shared-model and arbitrary-order boundaries, or a stricter M0 threshold for an explicitly identified exon, conditional on the same empirical interval family. It is not a first discovery of context effects, a calibrated measurement error, causal nucleotide interaction, marker-invariant intended-edit effect, full-haplotype assay, independent replication or inclusion/abundance separation. The increment is a complete source-specific comparison of the two exact compatibility boundaries between the already tested grid points. A null increment must be reported as such; no broader biological analysis is authorized.
