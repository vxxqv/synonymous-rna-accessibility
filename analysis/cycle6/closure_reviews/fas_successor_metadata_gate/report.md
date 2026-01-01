# FAS successor C: exact-vertex measurement metadata gate

**NO-GO for the proposed eight-vertex measurement-discrimination test using currently authenticated earlier FAS material.** Only five vertices, including WT, are in the declared 2016 analyzed single/double-mutant design. Its pooled assay measures included-RNA enrichment again. The separately validated clones are single mutants, whereas none of the seven non-WT certificate vertices is a single mutant. No genotype-matched inclusion/skipping record for those seven variants was recovered. This is a source-contract limit, not a biological negative or a claim that an undiscovered record cannot exist.

## Charter

Goal: decide whether existing earlier FAS measurements can separate inclusion from abundance/recovery for the fixed certificate vertices `{13,16,35,39,42,46,51,70}`. Scope/object/subject: ambiguity of included-RNA enrichment / a separately generated same-exon measurement check / exact sequence membership and source measurement/preparation records. This is a descriptive metadata gate. Only frozen sequence design and existing primary Methods were inspected. No 2016 outcome table, count, header, missingness or phenotype value was opened. This is a separately authored review with inherited context, not a blinded external review.

## Canonical sources

- Julien et al. 2016, [10.1038/ncomms11558](https://doi.org/10.1038/ncomms11558), [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4866304/), Methods: Doped library construction, Cell transfection, Validation of selected mutants, Computation of splicing enrichment scores, Data availability. Existing primary XML was reused; no new data retrieval was needed.
- [PRJEB13140](https://www.ebi.ac.uk/ena/browser/view/PRJEB13140), identified by the article as its sequencing accession. No raw/processed observation files were opened.
- Frozen cycle6 sequence contract SHA256 `f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38`. Vertex indices below are zero-based, matching the frozen certificate. Native reference mapping and intended included-protein identity were authenticated in the earlier sequence audit; actual minigene protein output remains unmeasured.

## Findings ledger

| Vertex | Substitutions from native human WT | Distance | Declared 2016 analyzed design |
|---:|---|---:|---|
|13|T26G, T35G, T39C, G41C, A44G|5|Outside radius2|
|16|T26G, T35G|2|Potential double mutant|
|35|T32C, A44G|2|Potential double mutant|
|39|T35G, T39C, A44G|3|Outside radius2|
|42|T39C, G41C|2|Potential double mutant|
|46|Native human WT|0|Reference|
|51|C24T, T26G, T32C, T35G, T39C, A44G|6|Outside radius2|
|70|C24T, T26G|2|Potential double mutant|

These are exact **theoretical design memberships**, not observed measurements. The paper's overall double-mutant completeness cannot be transferred to these four doubles without an authorized presence/coverage check. The initial doped pool can contain higher-order genotypes, but the published downstream analysis excludes variants with more than two substitutions. Their possible existence in raw material is not an authenticated processed-data contract and does not authorize raw-data expansion.

### Assay and preparation map

| Material | Source-supported unit and measurement | Fit to this successor |
|---|---|---|
|Pooled 2016 library|Triplicate HEK293 transfections. The exon-included RT-PCR band is excised and sequenced; variant frequencies are normalized to plasmid-input frequencies. Three independent input amplifications are technical preparations. Unique eight-base tags identify preparation libraries; two paired-end sequencing lanes are not extra biological replicates.|At most WT plus four theoretical double-mutant matches; actual presence remains unread. Same included-RNA proxy does not separate inclusion from abundance/recovery.|
|Individually tested clones|25 selected **single-mutant** minigenes, Sanger-verified, transfected in triplicate in HEK293 or HeLa. Inclusion and skipping gel bands are quantified to estimate their ratio.|No non-WT target vertex belongs to the documented single-mutant category. Their seven paired inclusion/skipping clone records are not recovered by this method description.|
|WT/control material|Native WT is the reference and basal inclusion is discussed in the source.|Exact clone/preparation-level paired-band record and uncertainties have not been authenticated here. Do not infer a complete eight-vertex dataset from a WT baseline.|

The pooled experiment does not supply an authenticated per-variant skipped-RNA denominator. If the variable exon is skipped, its genotype bases are absent from that molecule. A pooled skipped-band count would therefore need a linked genotype barcode or independent clone identity to assign the denominator to a variant. Source preparation tags alone do not provide that mapping. No such linked record was authenticated in the inspected Methods.

A direct included/(included+skipped) measurement could help cancel total transcript abundance for the same clone/preparation. It would still need correction or bounds for differential isoform RT/PCR/band recovery. Neither a second enrichment measurement nor a correlation with a few selected single-mutant validations resolves the current eight-vertex ambiguity.

## Scope filter

Excluded2016 outcome lookup, raw-read reconstruction, other datasets, new experiments and model fitting. No candidate was replaced after seeing current outcomes. Theoretical five-of-eight overlap was not relabeled observed coverage or sufficient independent validation.

## Object filter

Kept the exact eight current certificate vertices and the original measurement-discrimination goal. Did not substitute the five overlapping enrichment-proxy genotypes for the required paired inclusion measurement. Kept same-exon/source dependence explicit: 2016 precedes 2019, is cited by it and can inform its design/interpretation. Separately generated transfections do not create an independent exon or an untouched validation source.

## Conclusions

1. The current source Methods do not authenticate a full eight-vertex dataset that distinguishes inclusion from abundance/recovery. Seven non-WT vertices are all multi-mutants; the documented direct-PSI clone validations concern single mutants.
2. A 2016 same-proxy overlap analysis would test limited same-exon enrichment reproducibility, not the selected measurement question. It is not an acceptable automatic fallback or reason to open 2016 outcomes now.
3. Stop this successor at the metadata gate. Retain the completed current pilot's conditional ordinal contradiction and absence of an envelope-certified barrier, with the unresolved measurement ambiguity.

## One next-test specification and missing record

**Test C, conditional on recovery of an existing eligible record only:** compare direct inclusion/skipping measurements for all eight exact sequences against the current certificate's ordering, using matched WT controls and separately documented preparations. Do not use the original data as training or select a subset that agrees. The endpoint, error treatment, joint-order decision, exposure status and validation contract would need a separate freeze before any future outcomes.

The single most informative missing item is **one existing sequence-linked per-clone/per-preparation inclusion-and-skipping record for these eight sequences**, with: exact clone sequence/identifier; included and skipped band/read quantities or calibrated PSI generating records; culture/transfection/library identifiers and pairing; same-condition WT control; RT/PCR or band-recovery calibration; and sufficient replicate/error documentation. A table of mean enrichment alone is insufficient. The record has not been recovered, and this task does not assert it was deposited. No author contact, new source search or new experiment is proposed.

**Acceptance gate:** all eight sequences and corresponding independent preparation identities authenticated, with direct inclusion/skip or independently calibrated inclusion and a defensible error contract. **Stop rule:** a missing target, only the same enrichment proxy, untraceable clone/condition pairing or unsupported uncertainty means no executable biological test. It is a failed source/measurement gate, not evidence supporting or refuting the biological ordering.

`membership_receipt.json` and the short reproducible sequence-only checker record all exact sequences, memberships and hashes. `source_manifest.json` pins the Methods source and exposure boundary. No observed 2016 value or measurement presence was inspected.
