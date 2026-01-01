# FAS pinned count-file format and source-import audit

**Verdict: the 12 pinned count files are headerless at their first records. All 12 first records are the native human WT.** The original frozen parser's header stop was appropriate for its declared contract. A transparent format-only amendment can retain all data records and preserve the predeclared primary 27-scenario ratios. The literal author import cannot be called faithfully reproduced under that corrected parsing.

## Scope and exposure

Root reported that frozen revision `7d8faf59` stopped with `Unexpected count header` before analysis. This audit does not change that run, its source files, or its frozen code. Only first physical line format and sequence identity were semantically inspected. File sizes and full-file cryptographic hashes were checked for provenance; later records, coverage, genotype missingness and expression values were not inspected.

Root had already exposed the first WT count in BR1 and TR1. This audit did not display numeric count values from any file. The script classified the first line's three fields, verified the declared sequence pattern and ancestral mutation annotation, and recorded only format/identity booleans and hashes. Raw bytes were read for full-file SHA256/Git-blob authentication, with no later-row parsing. No statistical calculation, graph scoring or phenotype model ran.

The first audit-script attempt used the unrecognized codec spelling `utf8-sig`, failing before field interpretation. That script and failure receipt are retained. The corrected spelling `utf-8-sig` completed the same bounded audit. No source data or frozen original was modified.

## Exact first-record contract

`first_record_receipt.json` records every label BR1-BR9 and TR1-TR3, local file name, byte length, SHA256, verified source Git blob and first-line hash. For all 12:

- The first line is not `Genotype<TAB>Occurrences<TAB>MutationCount`.
- It has exactly three tab-separated fields matching declared sequence / nonnegative integer / nonnegative integer.
- Its sequence is the native human WT, not the inferred ancestral exon.
- Its mutation annotation matches distance from the source ancestral sequence.

This establishes the headerless first-record format; it does not claim a later-record or completeness scan. The source's intended collapsed-count contract has one row per genotype. Any later duplicate discovered by the amended reader must still stop under the frozen duplicate policy.

## Authenticated cause and literal author consequence

The upstream count generator `004_count_mutations.py` emits a textual header. However, pinned author `005_remove_sequencing_errors.sh`, blob `ba3858d75e4072c723c53c4842a9f88a792aa490`; the exact source content and raw retrieval are retained in the preceding source audit as `source_author_pin_1.txt` and `author_pin_1.json`. Its command filters with the 63nt genotype regex using `grep ... > ...3072.counts`. The textual header does not match that regex and is therefore removed. The first data row is retained by this filter.

Pinned `003_HEK293_Phylogeny_Library_Analysis/001_enrichment_scores.R`, blob `87838a729c4304ce3e38df87daaefb30c249119e`, then calls `read.table(This.File.Name, header=T)` for both BR and TR filtered files, before replacing column names with Sequence, Counts, Mutations. Under that literal call, the first physical line is consumed as column names; assigning new column names afterward does not restore its data. Thus each observed first WT record would be lost.

Under the intended one-row-per-genotype source contract, WT would then be absent in each imported table and absent from the TR1 row index. A WT-relative score cannot be computed from those literal imported tables. This is a conditional code consequence, not an assertion based on an unauthorized later-row duplicate scan. The intended source collapse supplies the uniqueness contract; the future parser must verify it.

This does **not** establish that the publication's released processed estimates are wrong. Original working files might have carried a restored header, import code might have differed, or later intermediate objects might have another provenance. Those possibilities are unverified. No processed ES, PSI, RData or publication outcome table was opened to decide them. Do not silently infer a historically successful run from the presence of analysis scripts.

## Minimal amendment recommendation

1. Preserve the original stopped freeze and receipt as failed source-format compatibility, without rewriting its status or hashes.
2. Create an explicit operational format amendment accepting the authenticated headerless layout. For these exact pinned files, the first WT line is data and must be retained. Prefer a strict format classification: accept either the exact known header or a valid first data record, record the mode, and reject an ambiguous/malformed first line. Do not make every arbitrary line skippable. Preserve all existing integer, mutation-annotation, duplicate, source-blob and genotype checks.
3. Keep the primary endpoint and all 27 separate-TR ratios unchanged. Retaining the source's first WT observation is required to implement that already declared endpoint; no threshold, model, cohort or favorable subset is changed. Do not infer anything from its magnitude. No coverage/value-based amendment is authorized.
4. Label the separate median-input calculation **“headerless-data reconstruction of the author's normalization formula”**, or equivalent. It is not the literal header-consuming historical import. Retain an explicit historical-reproduction gate: literal published code on these bytes drops the WT record, so a historical WT-relative point is unavailable under that contract. Do not fabricate a NaN/undefined source value as zero or use another reference.
5. Separate corrected all-record/TR1-aligned totals from the literal first-record-dropped totals if that historical diagnostic is ever authorized. No such total computation was performed here. The primary relative-ratio formula still cancels totals algebraically; it does not need a change to normalization.
6. Before any further real analysis, add synthetic tests for exact headered and exact headerless input equivalence, first-record retention, an invalid first record, duplicates, mutation annotation and source hash mismatch. Freeze the minimal adapter and revised labels/provenance before outcomes are parsed. Do not reopen source 2016 validation or retune anything after this source-format discovery.

## Decision boundary

**GO for a minimal format-only successor after review and freeze. STOP any claim of faithful historical WT-relative reproduction from the literal author import as currently documented.** The primary observational question remains the same conditional enrichment-order/72-state accessibility question. This audit adds no evidence about whether its null will be rejected, whether any functional obstruction exists, or whether such a result is biologically calibrated.

Source paths, immutable file hashes and read boundaries are in `audit_receipt.json` and `first_record_receipt.json`. This is a separately authored audit with inherited context, not a blinded external review. No external write, Git operation, source modification or new statistical analysis occurred.
