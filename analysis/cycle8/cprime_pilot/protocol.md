# B1: finite marker-context sensitivity and shared ordering

Status: pre-outcome specification awaiting independent implementation review and an immutable freeze. No source outcome may be decoded by this protocol until root records authorization against the final hashes.

## Question and biological object

For the 30 reference-annotated synonymous intended edits with at least two canonical marker states in Supplementary Data 8, does the finite observed marker set admit a common response within the empirical preparation ranges? How large can or must its multiplicative dispersion be? As a secondary context check, can one unrestricted latent ordering of those same local states explain both cell-line interval sets through separate nondecreasing response functions?

The object is exact local search-motif representation in RNA relative to DNA. It is not a full realized haplotype, total RNA, a protein amount, an isolated intended-edit effect, a clinical effect, or an edit-marker interaction. Source-described marker non-neutrality is established prior. The increment is a complete, uncertainty-explicit finite-set audit with exact bounds and all outcomes retained. No generic mathematical novelty is claimed.

## Fixed ascertainment and source

Use DOI10.1038/s41467-026-78344-9, Supplementary Data 8 sheet `Supplementary_Table_3_Multi_exo`, and the pinned author tree30125415304664c33370cd7801f820a345834e40. The existing B0 source contract governs annotation and sequence aliasing. Include `synonymous_variant` and `splice_region_variant&synonymous_variant` annotations. Use the 335 canonical IDs present in the metadata table, not all341 design aliases. Eligible groups are the 30 intended IDs with multiple canonical markers, comprising66 distinct local states. Preserve both cells, KBM7 and K562, and all three source rows per state/cell. The other33 single-marker synonymous groups are design-ineligible for a marker contrast, not outcome exclusions.

The arrayed R1207R panel is excluded from this pilot because its T-marker identity and matched WT endpoint are unresolved. That panel's discrepancy is preserved in B0; it is not repaired or used as validation. No author-predicted bad marker is removed from the multi-exon set. Selection is entirely metadata-based. IDs, motif strings, source hashes and expected row multiplicities must be sealed before values.

## Endpoint and source checks

Read only eligible IDs, cell/exon metadata, paired `RPM_RNA`, `RPM_DNA`, and `cprime_score` for source reproduction. Preserve each source row's RNA/DNA pairing. Do not infer pairing between different marker states or between cells from row order. Compute q=RPM_RNA/RPM_DNA without another pseudocount. Source code derives RPM from (count+1)/(total_reads+1); the processed table does not authenticate the raw count/denominator likelihood. The published score is a reproduction check, not an alternate outcome selected if q is inconvenient.

Treat deposited decimal strings as exact rationals for order and interval calculations. For the independently rounded published score, require absolute disagreement from q no greater than 1e-8 times max(1,abs(q)); this is a file-consistency tolerance, not a biological tolerance. Preserve the difference for every row. For positive paired RPM values, compare every finite source score, including a nonpositive score, before endpoint-usability filtering. An inconsistent score is a source-reproduction failure and prevents interpreting the pilot. Do not change the cutoff after inspecting failures.

Each state must have exactly three metadata rows in each cell. Missing, nonfinite or nonpositive RNA/DNA values or scores are recorded explicitly. Any missing required measurement value makes that state's primary sensitivity interval (0,infinity), retaining the state. A state with three valid rows has closed interval [min q,max q] and central median q. An absent or extra metadata row is an integrity failure because the sealed metadata require exactly three rows per state/cell; it is not recoded as ordinary missing measurement. No outcome-dependent state or group deletion. Retain all rows, missingness reasons, and complete/incomplete counts. Nonpositive values are not silently logged or treated as physical zero concentration.

## Exact finite-set quantities

For a group with M>=2 marker intervals I_m, define F(y)=max_m y_m/min_m y_m for strictly positive y_m. Primary lower bound is max(1,max_m L_m/min_m U_m), with unknown intervals interpreted as (0,infinity). Primary upper bound is max over DISTINCT m,n of U_m/L_n. Infinite bounds and unattained suprema are reported explicitly. The log2 width is display only; no floating logarithm decides a tie or a claim. The independent derivation must establish attainable/supremum conventions.

Report all60 intended-edit/cell combinations. The main descriptive counts are F_min>1, F_min>=2, F_max<=2, and unresolved/missing. The factor2 comparison is an operational display threshold, not a validated biological significance cutoff. Also report continuous exact bounds and the fixed1.25,1.5,2,4 fold grid so no selected threshold supplies the conclusion. Feasible common response means only that the interval intersection is nonempty. It does not establish marker neutrality or prediction of new markers.

Baselines: the point-median marker dispersion and the reference-assumption prediction F=1. The first can overstate certainty and the second is a testable simplifying assumption, not a fitted model. Do not report a p-value, calibrated confidence interval, population prevalence or predictive superiority.

## Secondary shared-order model

Within each of the five metadata-eligible exons (17,25,27,39,42), consider all eligible canonical states and the two cell contexts. The null permits an arbitrary real scalar z_v for each local state and separate nondecreasing functions h_cell,exon, with h_cell,exon(z_v) inside that state's empirical interval. This is a deliberately capable necessary condition for any shared sequence-additive model, not itself an additive fit. Cross-exon numerical orders are not compared because assay normalization differs.

For each cell, Ui<Lj forces z_i<z_j. Union the directed constraints from the two cells. An exact directed cycle is an infeasibility certificate. If the union is acyclic, a topological order and monotone interval completion in each context witness feasibility. Missing-state intervals contribute no unsupported order. Preserve every edge, source cell and rational inequality, every detected certificate, and feasible witnesses. Feasibility does not prove additive representability, physical mechanism or transport to unobserved variants. Infeasibility excludes a common scalar order only conditional on these descriptive envelopes and source identity; context-dependent normalization/recovery remains a competing explanation.

This secondary test is not called codon-additive validation. A full nested additive/context-specific/interaction ladder remains gated until its coding objects and design ranks are authenticated. No model fitted after seeing these outcomes will be confirmatory on these data.

## Dependence, sensitivity and validation scope

The three anonymous rows are an unordered set of source-reported preparations. Their range is a descriptive sensitivity set, not a confidence interval. Shared preparation or normalization offsets can change range geometry; pairing-free does not mean independence-free. A single common positive scaling within an exon/cell leaves within-context ratios and order constraints invariant. Preparation-specific offsets are not known and cannot be assumed to cancel.

Primary envelopes retain all three rows. A declared widening sensitivity uses [L/gamma,U*gamma] for gamma=1,1.25,1.5,2, with unknown intervals still (0,infinity). Recompute the same dispersion and shared-order tests at all four factors. This is a specified uncertainty family, not an estimated correction. No omitted-row or predicted-marker filtering search.

There is no fitted predictive model or data-dependent criterion in this audit, so an artificial split cannot provide validation. All groups are analyzed once after freezing. Cross-cell comparison is a second context from the same study and shared methods; it is not independent-study validation. The source's published examples are exposed prior evidence. Any successor prediction study must freeze a new question and whole-intended-edit holdout before its outcomes, and this audit will have exhausted these rows for fresh confirmatory reuse.

## Acceptance and stopping

Before values, verify metadata coverage and exact source hashes, synthetic positive/negative/tied/missing/zero cases, distinct-marker bounds, row permutation invariance, within-row pairing dependence, shared-order feasible witnesses and cycle certificates, and poisoned unselected cells. A fresh reviewer must inspect the implementation. The production reader requires the frozen protocol/code/input hashes and root's authorization and writes only a new output directory.

After one run, preserve all60 group/cell records, all132 state/cell records and all396 source rows, or their explicit missing-state replacements, the five exon comparisons, all widening sensitivities, every failure and the independent numerical reconstruction. No number of positive cases is required. Source-reproduction failure is a data-contract failure; wide intervals are information limits; feasible common responses are non-exclusion. No positive result triggers automatic broad expansion.
