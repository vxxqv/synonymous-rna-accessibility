# cPRIME selected-row provenance follow-up

**Verdict: the source lineage is partly recovered, but authoritative assignment of the 396 published rows to preparation labels and RPM totals remains unresolved.** No new phenotype, count, RPM, outcome-missingness or notebook-output values were inspected. The completed frozen analysis is unchanged.

## What is now authenticated

The pinned author revision is `30125415304664c33370cd7801f820a345834e40`. Cached notebook and counting-script bytes match the Git blob identifiers in its repository inventory. An additional ordinary public retrieval of `multi_exon_screen/analysis/gmm.ipynb` succeeded once and matches its pinned blob. Only code-source JSON fields were decoded; notebook outputs were skipped. No author analysis was executed.

The decisive source is [analyze.ipynb at the pinned revision](https://github.com/Boettcher-Lab/cPRIME/blob/30125415304664c33370cd7801f820a345834e40/multi_exon_screen/analysis/analyze.ipynb):

| Source location | Evidence | Consequence |
| --- | --- | --- |
| Notebook cell 6; cached source-only lines 42-68 | Reads `*_results_total.tsv`; derives modality, exon, cell line and replicate from filename tokens; computes `(counts_complete+1)/(total_reads+1)*1000000`. | The upstream computation has explicit sample-label keys and a defined processed-count denominator. |
| Cell 7; lines 71-78 | Pivots RNA and DNA using combined variant, exon, cell line and replicate. | Pairing is by recorded labels, not by positions of independently sorted RNA and DNA rows. |
| Cell 10; line 164 | Writes `wide_annotated_new.tsv` before dropping columns. | This labelled intermediate could close the row-identity gap if authenticated and linked to the published workbook. It is absent from the pinned repository. |
| Cell 10; lines 166-177 | Selects a 44-column `savedat` table without `replicate`, then writes `rpm_new.tsv`. | The ordered column list exactly equals the published `Supplementary_Table_3_Multi_exo` schema. This is the concrete label-loss point. |
| Both deposited analysis notebooks | No export to `Supplementary_Table_3_Multi_exo` or to Supplementary Data 8 is recorded. The second notebook reads a replicate-averaged intermediate. | The final TSV-to-workbook transformation and whether it preserved row order are unrecorded in the inspected source. |

The source-selected metadata comprises 396 worksheet rows, 132 state/cell triplets, 66 local sequence states, 30 intended edits, and exons 17, 25, 27, 39 and 42. These are identity counts, not usable-outcome counts. Identity-only extraction used columns A, B, E, F and G. No other cell payload was decoded.

Official cached [ENA metadata for PRJNA1257064](https://www.ebi.ac.uk/ena/browser/view/PRJNA1257064) contain 72 unique multi-exon library names covering both modalities, six exons, two cell lines and R1-R3. Sixty libraries, or 30 source-labelled RNA/DNA pairs, are relevant to the selected five exons. `labelled_library_pairs.json` gives their exact SRR and SAMN accessions. `selected_row_candidate_map.json` retains all 396 worksheet row numbers and the three admissible label pairs for each row. It deliberately sets the authoritative row assignment to null.

## Why this does not establish a row-level pairing

Within each state/cell triplet, all approved identity fields are identical. The source export removes the only R1/R2/R3 discriminator. A permutation of those three labels changes no published identity field. Consequently, matching the 44 headers and the repeated metadata does not distinguish the six possible assignments within a triplet.

The notebook supplies a plausible conditional reconstruction: if the exact `rpm_new.tsv` was inserted into the workbook without reordering and the original pivot/merge ordering is established, its retained row order may reproduce the replicate labels. That is stronger provenance than an arbitrary guess, but the missing intermediate and final export history prevent it from being an authenticated assignment. No row was labelled R1, R2 or R3 from its ordinal position, and no numerical fingerprint matching was attempted.

The RNA/DNA library-name pairing is itself a source naming convention. Distinct BioSample accessions do not prove distinct cultures or quantify shared-preparation covariance. Common exon and replicate names also do not make different exons independent biological units. The source states biological replication, but this pass does not reconstruct culture-level preparation identities.

## RPM totals are a separate missing record

The [counting script](https://github.com/Boettcher-Lab/cPRIME/blob/30125415304664c33370cd7801f820a345834e40/multi_exon_screen/count_edits_in_fastq/count_edits_fast.sh), lines 133-165, uses FLASH merged fragments and the quality-masked merged file for motif counting; the optional concatenation including unmerged reads is commented out. `count_ss.py`, lines 41-61 and 84-89, increments and exports `total_reads` for that processed stream. `merge_output.py`, lines 38-61, carries that total into the pseudocount RPM calculation. The count wrapper also writes merged/unmerged diagnostic totals in result-file comments, lines 214-220.

Those processed result files and their totals are absent from the pinned tree. Supplementary Data 8 does not retain totals or counts. ENA accession metadata cannot substitute for the number of quality-processed, FLASH-merged fragments; a raw archive read total would be a different denominator. No RPM total was inferred from a displayed ratio, raw-read inventory or assumed library size.

The upstream notebook also prints a warning rather than stopping if multiple result files match a directory. Its pivot aggregation and actual file inventory therefore need authentication before asserting one published RPM is exactly one archived run's unaggregated measurement.

## Exact closure objects and limits

An author-deposited `wide_annotated_new.tsv` preserving `replicate`, plus an authenticated row-preserving transformation to the published sheet, would resolve the selected row labels. A row-specific supplementary key mapping directly to named RNA/DNA libraries would also suffice. To recover RPM denominators, the required record is each labelled `*_results_total.tsv` or an equivalent processed-count manifest with its exact filtering/merging provenance and totals. Culture/preparation linkage requires explicit sample metadata beyond library names if it is to support covariance or independence claims.

These objects were not found in the pinned source tree or either analysis notebook. This bounded pass did not inspect another data cohort, read raw sequencing files, contact authors, repeat known access denials or change the frozen endpoint. The existing unordered, permutation-insensitive interpretation remains the supported scope. A paired covariance or library-depth extension remains gated.

## Receipts and exposure

- `lineage_receipt.json`: source blob checks, exact identity counts, schema comparison and output hashes.
- `selected_row_candidate_map.json`: all selected worksheet rows, with unresolved assignments explicit.
- `labelled_library_pairs.json`: all 30 source-labelled candidate pairs and their official accessions.
- `gmm_retrieval.json`: URL, pinned blob, SHA256 and source-only parsing boundary for the sole new retrieval.
- `audit_lineage.py` and `trace_notebook_source.py`: inspectable metadata-only checks.

No current frozen numerical result was opened. An existing Methods text search incidentally displayed published GMM summaries already contained in the cached article, including the BIC and classification-concordance prose. These were not used for this provenance decision or treated as new analysis. No new variant-level values or missingness were examined. Full papers, notebooks and workbook evidence caches are not automatically approved for redistribution.
