# Synthetic preflight: FAS count reader and normalization

**Final verdict: approve the corrected numerical count-reader contract narrowly, with the nonblocking output-provenance recommendation below.** The initial snapshot required revision: one valid missing-data case incorrectly aborted the primary analysis, and nucleotide/codon equivalence needed explicit interpretation. Both were corrected and checked before observed count access. Neither finding uses observed count values or source coverage. Historical failing snapshots and results are retained.

The held `count_analysis.py` snapshot has SHA256 `b0ee9e7282f6cffae5f48be1739a04167093c093650eb4fe94a8afac7cb9e7ce`. The reviewed draft is `731620b9d8d4057053b80143b2c476b94f8c3feb46c6342c9841ff531ad3393d`. Immutable copies accompany this review. Shared source was not edited. Review is separately authored with inherited context, not blinded external review.

## Required corrections

### R1. A zero TR1-aligned total must not erase a valid partial primary analysis

Location: initial `source_relative`, line49, called by `analyze` after constructing the primary envelopes. The unconditional `if any(t<=0 for t in aligned): raise ValueError(...)` treats a missing/undefined author normalization as an invalid entire primary dataset. That conflicts with the fixed full-domain unknown-state contract.

Minimal synthetic input uses only valid declared genotypes and nonempty files. Let v be any selected genotype other than WT, and a be the declared ancestral genotype outside the native-reference class:

```
TR1 = {WT:1, v:1}
TR2 = {a:1}
TR3 = {WT:1, v:1}
BR1..BR9 = {WT:2, v:3}
```

TR2 is nonempty and has no malformed count or duplicate. Its selected/reference observations are missing, so the primary27-scenario envelope must mark affected vertices unknown. It is not permissible to silently remove TR2, but it is also unnecessary to abort all primary geometry. The author source aligns TR2 to TR1, yielding NA rows and a zero `na.rm` total. Its undefined input frequencies are omitted by `median(...,na.rm=T)`. The remaining inputs give source-relative WT=1 and v=3/2 for every output. The current code instead raises `ValueError('Empty author-aligned source universe')`.

Required fix: handle an aligned total<=0 as undefined normalized frequencies, preserving the source's NA/NaN omission semantics for input medians and recording the degenerate total. For an output with zero aligned total, preserve undefined source estimates for that output. Do not manufacture finite normalization or discard a primary scenario. Primary all72 envelopes, unknown diagnostics and any valid original-source reproduction must remain reportable separately. Add tests with no shared TR1 support and with explicit all-zero aligned values; the latter produces NaN rather than NA but `na.rm` still removes it in the source median. If every normalized input for a genotype is undefined, its source score remains unknown.

### R2. The nucleotide sensitivity is equivalent to the codon null here

Location: initial `analyze` feature construction, lines69-73; draft coverage/model-check paragraph. Independent exact SymPy ranks are:

| Matrix | Rank |
|---|---:|
| Eight-column codon categorical design |8|
| Eight-column nucleotide categorical design |8|
| Concatenation of both designs |8|

Each Leu factor occupies three of the four combinations at its two variable nucleotide coordinates. Intercept plus the two binary nucleotide indicators spans all functions on those three states. Combining the two Leu factors and three binary codons preserves equality of the complete latent column spaces. This equality is domain-specific; it is not a general statement for larger genetic-code fibres.

Therefore the nucleotide model is not a weaker comparator on this fixed graph, and a P3 middle-state dip can arise even from nucleotide additivity restricted to that missing-corner domain. Keep the code only as a redundant exact-space/solver consistency check, or label the sensitivity's equivalence explicitly. Any disagreement in certified feasibility would indicate numerical/certificate handling or implementation inconsistency, rather than evidence for within-codon effects. Do not describe a codon-versus-nucleotide contrast as a biological result here.

## Checks that passed

`independent_checks.py` uses held snapshots and independent Fraction arithmetic with20 synthetic source fixtures. It does not invoke the author's functions as its numerical reference. The source-normalization reference explicitly forms normalized enrichments and their WT ratio. Cases include missing rows in a nonreference input, missing output rows, and genotype representation outside the TR1 index.

-12960 source output-by-vertex values agreed exactly where the source normalization was defined.
-38880 separate-input relative ratios agreed exactly.
-72 zero-output-reference cases and72 zero-input-reference cases correctly returned unknown; a finite variant-output zero correctly returned zero.
-Native human WT is used for relative scores; mutation-number validation correctly uses the ancestral source annotation. These are different reference roles.
-TR1 row alignment and `na.rm` omission are reproduced for ordinary missing rows. Source median of normalized inputs is not replaced by a median of reference ratios.
-The primary envelope becomes `[0,infinity]` when any of a vertex's27 scenarios is unknown. The genotype list is not reduced to observed rows.
-The author-point expression is the arithmetic mean of nine per-output relative scores, conditional on all nine being defined. It is not a pooled-count ratio or ratio of means. Its label `mean_author_per_output_relative_ES` makes that derived-point convention explicit; it is not claimed identical to the original published mean PSI normalization.
-The code checks that a finite author-relative point lies inside a fully finite scenario envelope. The general ratio-of-medians enclosure proof is in the preceding measurement review; finite checks alone do not prove it.

The supplied author test suite's constant fixture does not exercise the zero-aligned-total case, because all its selected rows and aligned totals are positive. Its passing1944 ratios and five guard cases therefore do not close R1.

## Nonblocking execution/provenance recommendation

The CLI writes the output path with `write_text` and does not refuse an existing destination. Use an explicit new-output guard before any scientific reads/output creation, or demonstrate the frozen orchestrator enforces the same requirement. This protects historical failure receipts; it does not change numerical results. Freeze exact source/implementation/protocol hashes after corrections. No actual source count should be opened merely to decide whether R1 occurs.

## Scope and stopping

This preflight covers parsing semantics, references, normalization, unknown propagation and feature-space interpretation. It does not replace the separate exact ordinal/bottleneck proof audit or numerical solver review. No new model, threshold or source was selected. The calibration literals remain masked in the preceding review and were not opened here. No count file, phenotype value, outcome missingness, historical2016 validation observation or raw read was accessed.

### Correction verification

The corrected count-reader snapshot SHA256 is `307c3a3fe27720bceeffb4a348d9b10b7cc048e2a47a19049ab529a23454edd4`; corrected protocol SHA256 is `eb6e76ed299f27227f68c541b9a0a3af2aee1792ca3efac239e1a46fe211eb51`. The source-normalization function now omits undefined input frequencies, records undefined aligned-total indices and returns unknown output estimates when required. The disclosed minimal example now matches the independent reference exactly, without an exception. The code verifies exact design-span equivalence and requires matching certified feasibility status for the two encodings; the protocol labels that equivalence rather than a biological sensitivity.

`independent_checks_v2.py` repeats the same fixed20 synthetic cases against corrected held snapshots: all12960 source comparisons,38880 relative comparisons and reference-zero checks pass, with no discrepancies. It adds an explicit mismatch failure in the formerly unreachable successful-degenerate branch. Only snapshot/output names and that disclosed verifier guard changed; the earlier failing verifier and receipt were not rewritten. This bounded delta review does not claim another independent proof or real-data replay.

R1 and R2 are closed for these exact corrected snapshots. Preserve the failing snapshot and receipt as historical evidence. One frozen descriptive pilot remains reasonable under the earlier conditional assay-measurement limitations; the CLI output-provenance recommendation remains for the final runtime freeze or orchestrator.
