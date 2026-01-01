# Shanas2026 metadata/code contract pilot

## Decision

**NO-GO for the requested authenticated same-protein/synonymous accessibility extension.** The author deposits authenticate the same **cassette allele index** across intended promoter contexts, not matched full coding sequences or invariant isoform proteins. All 9,608 annotated cassette-exon lengths are 2 modulo 3. The reference used for read mapping is a partial reporter amplicon; unused legacy full-reporter strings do not establish the 2026 coding start, frame or terminating stop. Synonymous edge/square counts are **not identifiable and were not computed**, not zero.

**Conditional metadata feasibility only** remains for a separately defined, single-endpoint matched-cassette RNA-context portability study. This would change the object from synonymous protein-preserving routes to context transfer of a measured RNA profile. The original paper already establishes promoter effects, sequence-dependent sensitivity, cryptic isoforms and predictive features. Family/context-heldout calibration is a possible methodological difference, not an established substantial biological novelty. Do not open phenotype data merely to see whether it looks promising.

Root should close this bounded test now. A next phenotype analysis is not recommended under the current same-protein contract. Any further action requires an explicit scope decision and a pre-outcome contract resolving the coding/isoform, preparation and endpoint issues below.

## What was inspected

Primary [Shanas et al.2026](https://doi.org/10.1038/s41467-026-76370-1) Methods/Discussion and the already retrieved ArrayExpress E-MTAB-16553 IDF/SDRF were reused. Both public author repository inventories were retrieved read-only. The mapping scripts and figure notebook's **source code only** were inspected. Author code was never executed.

Before any mixed table was opened, `selective_readers.py` specified the permitted fields: `library index`, `subset`, `gene name`, `exonlength`, `exonstart in sequence`, `exonend in sequence`, `sequence`. It discards unselected CSV fields while scanning and skips unselected spreadsheet cell values; phenotype values are not converted, retained, printed or used for missingness. Notebook outputs are skipped lexically without JSON-decoding them. The full notebook and mixed raw tables were not saved. Selected metadata and exact transfer/source hashes are retained.

The code-only notebook has 248 cells, including 230 code cells. Its upstream Git blob is `2100269266db24477ca212ba516ae957276a5e95`; transferred bytes SHA256 is `f139c91791ac524cffcac558abfee7c8d8563691abd01ea1ed0a316a4bfec425`. The retained source-only file has SHA256 `d58dca5d2e977693f705881362f8388a626a7d9e73404e3be8aa5f1686bb0021`. No output cell or deposited assay phenotype/missingness was inspected. Synthetic tests exercise notebook output rejection, CSV quoting and spreadsheet selection with an intentionally invalid unselected numeric field.

## Authenticated sequence/design inventory

| Metadata quantity | Count / result |
|---|---:|
| Unique library IDs | 9,608 |
| Sequence length | 210 nt for every row; A/C/G/T only |
| Unique full 210-nt strings | 9,608 |
| Unique 12-nt barcodes at `[18:30]` | 9,608 |
| Unique 147-nt variable strings at `[45:192]` | 8,924 |
| Nonempty gene labels | 134 |
| Rows without a gene label | 773, all in `cassettecomb` |
| `cassette_filtered` rows | 264 |
| Distinct gene labels / variable sequences in that subset | 132 / 132 |
| Records per such reference gene | Exactly two, with different barcodes |
| Annotated cassette exon length | 2–89 nt; all lengths are 2 modulo 3 |
| Annotated exon end in 210-nt source | 162 for every row |

The figure-repository workbook and mapper-repository CSV agree **exactly on every ID, sequence, subset and exon-coordinate field**. The only selected-field difference is representation of 773 empty gene names as empty CSV strings versus empty workbook cells; after empty-value normalization the gene labels also agree. The source-identity chain is therefore materially stronger than nominally matching filenames.

The paper describes an original 2×134 endogenous design. The current subset named `cassette_filtered` contains 2×132 records. This is a precise subset-count distinction, not evidence that four phenotype measurements were lost or that the publication is erroneous. No outcome filtering, coverage or missingness was read to explain it. A later cohort must use an authenticated metadata rule rather than assume 268 complete reference measurements.

The shared 147-nt strings and duplicated reference barcodes must remain grouped in any predictive split. Gene labels provide family grouping for 8,835 rows in principle, but 773 composite rows have no gene label; they cannot simply become 773 independent families. Sequence equality is known, not sufficient phylogenetic or biological independence. No actual folds or graph were selected here.

## Coordinates, intended RNA isoforms and coding limits

The used mapper reconstructs

`ref = sequence[:30] + const1 + sequence[45:] + const2`,

where the retained constants have lengths 288 and 833 nt. The resulting reference is 1,316 nt. Source sequence positions after 45 map with offset +273; the splice-detection script independently adds that offset to exon annotations. The notebook slices the same annotations as zero-based half-open intervals.

The paper's approximately 800-nt downstream constant and the mapper's stale `1839` comment must not replace the actual 833-nt string. Exact source arithmetic settles the amplicon mapping while leaving the full-ORF question unresolved.

| Part | Mapping-reference interval / boundary |
|---|---|
| Barcode | `[18:30]`, retained upstream of the first constant splice donor |
| Variable region | `[318:465]` |
| Upstream constant donor | 88, followed by GT |
| Annotated cassette exon | `[exonstart+273:435]` |
| Cassette donor | 435 |
| Downstream constant acceptor | 1214, preceded by AG |
| Canonical included amplicon | `ref[:88] + ref[exonstart+273:435] + ref[1214:]` |
| Canonical skipped amplicon | `ref[:88] + ref[1214:]` |

These are intended RNA/amplicon coordinates, not full transcript or CDS annotations. The source code also recognizes alternative splice sites and retained introns; actual isoform frequencies are phenotype data and remain unopened. The presence of two splice junctions alone does not authenticate the native intended protein product.

The primary Methods place the library integration site at the GFP coding terminus **before the stop codon**. The used reference omits the complete upstream GFP coding prefix/start. A `const0` string and a longer `tot_seq` string exist in the source, but are not used to construct the mapped reference. The splice-detector source retains an older extended tail, while the mapper has a shorter edited tail. Those literals are not a versioned, annotated 2026 full-plasmid/isoform reference for every promoter construct. Taking an unused literal as the actual full reporter would be an unsupported substitution.

The modulus result is decisive against a naive coding interpretation: if the cassette is translated continuously, inserting/skipping an exon of length 2 modulo 3 changes downstream phase. If a stop upstream renders it untranslated, that must be authenticated instead. The current evidence does not choose between those coding/termination possibilities or establish a fixed main-ORF product. Barcodes are also retained in the upstream amplicon exon, so ignoring them for protein identity needs a CDS-location proof. Translating an isolated cassette exon or the 147-nt segment cannot supply that proof.

Therefore even identical cassette sequences across promoter contexts do not prove identical actual mRNA isoform proteins. Identical intended full-CDS sequence, invariant translated allele product and invariant measured isoform mixture are three different assertions. None is replaced by matching the allele index.

## Preparation and context unit contract

The already authenticated SDRF has 100 FASTQ file rows corresponding to **50 source samples**: 40 targeted reporter and ten bulk RNA-seq sources. The doubled file rows are paired-end mates. They are not 100 biological preparations.

For the main five-promoter comparison, source names `pm425a/b`, `pm454a/b`, `pm455a/b`, `pm456a/b` and `pm457a/b` map to pEF1α, full CMV, short CMV, ACTB and minimal CMV, with labeled replicates 1/2. The main notebook dictionary and SDRF support that mapping. A second notebook dictionary (`samplenamedict_prom`) inconsistently labels pm456 as short CMV; the SDRF/main dictionary, not that isolated label, must govern any future contract. Treat it as a source-code labeling issue, not evidence of swapped biological samples.

The primary Methods explicitly describe duplicate cloning of each library and separate per-sample transfection/harvest procedures. The paper calls the repeated measurements biological replicates; the SDRF gives distinct sample/extract/run identifiers. This supports two labeled preparation chains per main promoter, but does not provide a complete clone-to-transfection-to-culture lineage or establish independence of all sources, shared stocks, batch effects and cross-context covariance. Duplicate cloning, read mates and replicate labels must not be multiplied into additional biological units.

Additional source names cover drug/control samples, upstream-intron/promoter variants and integrated/transient material. These are different contexts or batches, not extra interchangeable replicates. The `_R` suffix in the notebook is generated from a file-family rule for `ceadd` files; it is not, by itself, proof of an independent validation experiment. Integrated material changes construct backbone, genomic placement and selection history, so it cannot be treated as an exchangeable repeat of transient expression.

The matching gene library is intended to be assayed in each context. **Actual measured variant/context intersection is unknown here**, because extracting IDs from retained outcome tables would reveal outcome availability and was explicitly prohibited. The design upper bound is not an eligible measured denominator. No PSI, expression, read count, outcome missingness or target availability was accessed.

## Endpoint and marginal novelty

The primary PSI method counts exon-included reads over all passing reads after isoform/variant filtering; it is not automatically the canonical included/(included+skipped) ratio and is not a splicing rate. A future canonical-isoform fidelity endpoint requires exact junction and denominator definitions. This source does not authenticate genotype-aligned GFP protein measurements for every cassette allele: bulk FACS expression checks and bulk RNA-seq cannot be silently converted to a second per-variant function.

The original source-only classifier calls use `StratifiedKFold(n_splits=5, shuffle=True, random_state=42)` and `cross_val_predict` with no gene-family or context grouping at the inspected calls. Its labels are clusters derived from the multi-context profiles. A prospective family/context transfer test would therefore ask a different question, but must freeze a training-defined response/label rule rather than derive target clusters using every assessment outcome. This difference is a potential methodological extension, not proof of novelty or an accusation that the original analysis answered its own stated question incorrectly.

Known promoter effects, cryptic isoforms, exon-length/motif associations and feature-based cluster prediction are prior findings. A same-cassette, held-out-family/context prediction beyond an appropriate global monotone context map and simple composition/splice baselines could be distinct in principle. It would still be an author-library reporter comparison, not independent native-gene replication, synonymous accessibility, an RNA-structure mechanism or practical-impact evidence. No thresholds, predictor or split were selected in this pilot.

## Cheapest next decision and stop rule

1. **Do not open phenotypes under the same-protein accessibility claim.** Synonymous counts remain null/“not computed” in `evidence.json`.
2. If root wishes to preserve strict synonymy, the specific missing source is an authenticated annotated 2026 reporter construct and intended isoform-CDS map, including the GFP junction, barcode location, start and terminating stop. The explicitly identified promoter-sequence supplement was not recovered here because its public CDN host failed DNS resolution, not because access was denied or the file is absent. No repeated recovery attempts were made. Such a reference could revise the coding verdict; absent it, stop the branch.
3. Alternatively, root may explicitly choose a **matched-cassette RNA-context** scope. Before any values, that requires one precise endpoint, family grouping for duplicate/composite designs, a sample-lineage/assessment-material contract, and a comparison to the original cluster analysis. If independent assessment or a distinct question cannot be specified, close this route rather than consume its outcomes speculatively.
4. The completed Mordstein follow-up remains an unresolved reference-index/pairing gate. Historical summaries should say “follow-up completed; gate remains unresolved,” not suggest an unperformed generic retry or imply the unavailable source has been recovered. No new Mordstein retrieval was attempted in this task.

This is a source/contract limitation, not a biological negative. No new experimental design, sequence intervention, phenotype analysis, model fit, upload, author contact, commit or publication was performed. Only the assigned directory was written. Numeric outcome verification is not claimed.
