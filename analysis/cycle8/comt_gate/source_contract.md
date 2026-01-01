# COMT A0 source contract

Status: sequence and identifier gate passes; uncertainty and cross-assay preparation pairing remain unresolved. No variant phenotype, count, confidence-limit or FDR value was inspected. This is a source audit, not an analysis result or a frozen pilot.

## Authenticated sources

The final article is Hoskins, Rao, Tante and Cenik (2024), *Integrated multiplexed assays of variant effect reveal determinants of catechol-O-methyltransferase gene expression*, DOI [10.1038/s44320-024-00018-9](https://doi.org/10.1038/s44320-024-00018-9), PMID 38355921, [PMC11066095](https://pmc.ncbi.nlm.nih.gov/articles/PMC11066095/). The saved Europe PMC XML is the final article. The primary record links the earlier preprint as an update, not an erratum. Targeted DOI/correction search and the retrieved record found no separate correction; this is not a guarantee that none exists.

The article explicitly links [GSE246139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246139), [author code](https://github.com/ijhoskins/integrated_MAVEs), MaveDB experiment set urn:mavedb:00000661 and [Zenodo 10403931](https://doi.org/10.5281/zenodo.10403931). Author Git tree 7ee507ae0cef0dad1b42b06883bb8fe3116334fc is pinned locally. Zenodo v1.0.2 was deposited on 18 December 2023 and predates the final article. Its 57,934,893-byte archive matches published MD5 e9e177a991dbb972cec247fa717ba57c and local SHA256 e3fc127f268c10efbb1d1c736423c51bc5aa900e9e363d66cc8683490516cbe2. A final-paper label must not erase this version history. Every retrieval URL and hash is in retrieval_log.jsonl; local file hashes are in source_manifest.json.

The four public MaveDB API records resolve successfully despite the website search result reporting a missing experiment. They are a-1 mid polysome, b-1 heavy polysome, c-1 RNA and d-1 protein. The records link the preprint and have later metadata modification dates. No outcome-level equality between these current exports, the final supplementary tables and a regenerated author pipeline is claimed.

## Sequence and nucleotide identity

The complete 6,232 nt author FASTA equals the sequence in the authored SnapGene plasmid exactly after case normalization. Its expressed fusion coding region is positions 4653 through 6203, followed by a stop at 6204 through 6206, using one-based inclusive coordinates. Translation gives 517 amino acids without an internal stop. The construct contains the native exon-2 5' UTR, N-terminal Flag, COMT, C-terminal moxGFP and an IRES-mCherry element. A plasmid reference is not an authenticated complete mature transcript isoform inventory.

The MaveDB target sequence is the same 357 nt slice, plasmid positions 4794 through 5150, in all four score sets. This covers native COMT codons 40 through 158, including the unmutagenized gap. The actual target libraries are codons 40 through 74 and 136 through 158. The author removes eight amino-acid positions from fusion annotations to obtain native COMT numbering. All retained synonymous calls were independently applied to the full fusion and translated; none changes that full fusion protein. This verifies synonymy relative to the authored construct, not every endogenous COMT isoform.

One source annotation defect is retained: GFF transcript/exon end 6233 lies one base beyond the 6232 nt FASTA and SnapGene sequence. CDS, start, stop and target intervals are in bounds and consistent. The original files were not edited.

Join by exact target reference hash plus hgvs_nt, or by canonical plasmid nucleotide edits. Do not join by hgvs_pro: synonymous nucleotide variants collapse at protein level. The assay uses direct targeted variant sequencing, not a separate barcode-to-variant lookup. Sample/sequence and preparation identities remain essential. Multi-nucleotide codon replacements are distinct from SNVs.

## Header and identifier coverage

| Export | All variant identifiers | Synonymous identifiers | Synonymous codon sites |
|---|---:|---:|---:|
| RNA | 3463 | 128 | 54 |
| Mid polysome | 3317 | 128 | 54 |
| Heavy polysome | 3117 | 127 | 54 |
| Low-fluorescence protein proxy | 3050 | 126 | 53 |

There are 125 synonymous identifiers shared across all four exports at 53 codon sites: 77 at 33 ROI1 sites and 48 at 20 ROI2 sites. Of these, 113 are SNVs, at 52 sites. RNA and protein alone share 126 synonymous identifiers at 53 sites. The theoretical target codons have 133 distinct synonymous nonreference codon replacements across 56 sites. These are identifier intersections, not verified finite-outcome or quality-filter coverage. No zero/missing phenotype test has been performed. All HGVS identifiers parsed; there are no duplicate nucleotide identifiers within an export.

GEO has 157 samples: 156 named COMT libraries and one conventional ribosome-profiling sample. The archive has 158 variant-summary files: the 156 GEO COMT libraries plus two simulated random-forest training summaries, COMT_108_train and COMT_109_train. Headers expose POS, REF, ALT, codon/amino-acid annotations, DP, CAO and quality fields. Only headers were read. MaveDB exports have accession, hgvs_nt, hgvs_splice, hgvs_pro, score, score_ci_lower, score_ci_higher and score_FDR. Only the first four identifier fields were decoded.

## Measurement and normalization contract

RNA is total-RNA variant representation relative to total gDNA representation. Mid and heavy load are polysome metafractions 3 and 4 relative to total RNA. The protein proxy is enrichment in low-moxGFP relative-to-mCherry population 3 versus total gDNA. Positive low-bin enrichment is not positive protein abundance. The paper reverses the effect sign for a figure, but that does not authorize assuming every deposited export has already been sign reversed. A source-formula check must pin orientation before evaluation.

The assay does not measure absolute RNA molecules, protein molecules, translation flux or a calibrated mean fluorescence. Polysome occupancy can change through initiation, elongation or other processes. Cells for RNA/polysome readouts had 24 h induction, whereas cytometry used 21 h. These assays share a construct and cell system but are not identical contemporaneous physical measurements. Other cytometry bins exist, but neither an invented mean protein amount nor stochastic ordering is justified without authenticated gate boundaries and bin sampling fractions.

Author methods specify four independent stable cell lines per ROI, with some repeated growth, preparation and sequencing. Technical replicates are summarized before testing. They are not extra independent biological units. Within a sample and nucleotide position, the R code derives WT_CAO as DP minus summed variant counts, and normalizes (CAO+0.5)/(WT_CAO+0.5). It subtracts the median negative-control normalized frequency by nucleotide variant and removes nonpositive residuals. Shared denominators, count composition and negative controls induce dependence. An independent error bar per variant would discard that dependence.

The pipeline then uses ComBat, with assay model covariates for RNA/gDNA but not polysome/cytometry; ROI2 RNA/gDNA has no batch correction. Complete-case filtering precedes this step. Its empirical normalization uses information across variants. Exported scores therefore cannot automatically be relabeled as outputs from a training-only prospective preprocessing scheme. A fixed published-score prediction question can use them, provided its transductive preprocessing history is explicit.

ALDEx2 uses paired=FALSE, iqlr denominator and 128 Monte Carlo samples. Its standardized median group effect and marginal 75% intervals are not physical fold changes, biological replicate standard errors or a simultaneous joint uncertainty set. The published intervals do not supply cross-assay covariance or familywise reversal certification.

## Source exclusions and actual ambiguity

The source rules include target-region restriction, amber-only nonsense calls, log10 variant frequency above -5.8, a simulated-error classifier for SNPs, negative-control subtraction, complete observations across the relevant biological replicates, and median natural-log frequency above -10. The code explicitly drops samples 120, 130, 146, 145, 47, 57, 173 and 174 plus initial five-fraction polysome samples 7 through 36. Source diagnostics label additional outliers, but a diagnostic comment is not itself an exclusion. No new subset has been chosen here.

**Preparation crosswalk blocker:** 20 of the 47 FACS sample labels disagree between GEO and code. GEO calls COMT_42 through COMT_51 biorep2, and COMT_52 through COMT_61 biorep1; code calls them FACS_bio1 and FACS_bio2 respectively. Later repeat samples make it unsafe to infer a universal harmless swap. GEO bins P3/P4/P5/P6 are also relabeled code quadrants 1/2/3/4. This may be deliberate relabeling, but no recovered authoritative record resolves the stable-cell-line correspondence across assays or repeats. Exact rows are preserved in replicate_label_crosswalk.json. No source was silently preferred to create paired observations.

The missing record for paired statistical uncertainty is an authoritative sample-to-stable-line crosswalk reconciling these 20 assignments, growth/preparation repeats, assay identity and source bin names. The missing computational object for joint published-score inference is a joint resampling/normalization construction or draws that retain shared WT, negative-control and preprocessing uncertainty. Their absence blocks a covariance-calibrated assay substitution or variant-reversal claim from the current score exports alone. It does not erase the measured information or prevent a descriptive proxy-level pilot with explicit dependence bounds.

## Exposure

Published abstract, methods and source-code comments were read, including the already published statement of directional/divergent effects, a published coverage comment, and two hard-coded fluorescence reference constants encountered in code search. A source figure annotation named an abundance hotspot. These are prior/public exposures, not new variant-score evaluation. MaveDB metadata responses also contain an automatic aggregate scoreDistribution field; it was downloaded but not intentionally inspected or summarized. Score CSV numeric tails and archive count bodies were not decoded or analyzed. The accession/variant-presence intersection is now exposed and cannot later be described as an outcome-blind assessment of source coverage.
