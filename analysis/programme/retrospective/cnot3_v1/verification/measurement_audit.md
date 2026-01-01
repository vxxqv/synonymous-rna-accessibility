# CNOT3 measurement and processed-data eligibility

**Firm verdict:** proceed only to a frozen **conditional descriptive gene-level** analysis of deposited point measurements. Comparable control/CNOT3 half-life fields and measured CNOT3-association fields are publicly available. The data do **not** currently identify a calibrated absolute molecular degradation-rate contrast, a causal CNOT3 pathway contribution, biological-replicate uncertainty, or expressed-transcript-specific codon composition. Those inference gates remain closed. Their failure is not a biological null result.

The original processed24-column conversion-ratio matrix is also available and traceable to the time-course libraries. It could support a separately specified descriptive trajectory-model audit. Refitting it is unnecessary for the smallest point-estimate analysis and would require its own model, normalization, background and growth-dilution protocol. It must not be called reconstruction of the published half-life fit without those missing settings.

No gene-level numeric half-life, IP-association, conversion-ratio or quality outcome was inspected. Data-bearing sheets were restricted to row1 string headers and, after explicit root authorization, only the whitelisted GeneID column. Complete explicitly labeled description/legend sheets were read as measurement documentation under separate root permission. No sequences were acquired, no model was fitted and no validation subset was selected.

## Authenticated sources and acquisition correction

The official Entrez package was used for `esearch GSE268325[ACCN]`, `esummary` and `efetch` with db=gds. GSE268325 is a SuperSeries containing GSE268322(SLAM-seq) and GSE268324(selective ribosome profiling). Official public GEO FTP-over-HTTPS SOFT files returned200, without credentials, challenge workarounds or browser retries. The three SOFTs contained metadata and no sample-table rows. Saved raw receipts and metadata preserve the original source encoding.

| Source | Authenticated available object | Scope inspected |
| --- | --- | --- |
| [GSE268322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268322) | [Processed conversion workbook](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268322/suppl/GSE268322_CNOT3_SLAM-seq_to_C_convertion_ratio.xlsx) | `CNOT3_SLAM-seq` row1 headers; Gene_ID columnB only; full `Table Legend` documentation sheet |
| [GSE268324](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268324) | Official metadata lists `GSE268324_CNOT3_selective_ribosome_profiling_All_footprints.xlsx` | Series/sample metadata and listed filename only; this large processed footprint table was not downloaded or inspected |
| [Original study](https://doi.org/10.1126/science.adq8587) | [PMC11583848](https://pmc.ncbi.nlm.nih.gov/articles/PMC11583848/) | Measurement-method semantics supplied by the separate primary-prose audit |
| [Reanalysis](https://doi.org/10.1002/2211-5463.70277) | Correct [gene-level TableS2, s002.xlsx](https://pmc-oa-opendata.s3.amazonaws.com/PMC13398649.1/FEB4-9999-0-s002.xlsx?md5=e44bf579167ec8a242b182c0fd926faa) | `TableS2` row1 headers; GENE_ID columnA only; full `description` sheet |

The initially supplied cloud object `s003.xlsx`, MD5`187fef9b96c9fc60cf84a3b82964c946`, is authentic but is **TableS1 codon-pair data**, with sheets`PTMD` and`description`. Its headers include P-site/A-site codons and enrichment, not gene half-lives. The primary XML explicitly maps supplementary item0002/TableS1 to s003 and item0003/TableS2 to s002. The corrected s002 MD5 is`e44bf579167ec8a242b182c0fd926faa`, matching the retrieved bytes. The wrong-table receipt/schema remain preserved; it was not a corrupt file or an access failure. Original-study TableS2 is a sequence/guide/primer inventory according to the prose lane and is not treated as fitted half-lives.

## Exact sample and library map

All SLAM samples are labelled HEK293T. The source calls the conditions`sgNT`/control and`sgCNOT3`/CNOT3-knockout or depleted. GEO names guides sgNT1,sgNT2 and sgCNOT3-1,sgCNOT3-2 collectively, but does not assign a particular guide to a GSM/replicate label. A label`rep1` is not proof of guide1 or an independent biological culture.

| Chase time(h, explicitly documented) | sgNT rep1 / rep2 GSM | sgCNOT3 rep1 / rep2 GSM | Processed SLAM labels, NT then CNOT3 |
| --- | --- | --- | --- |
| 0 | 8289652 / 8289653 | 8289640 / 8289641 | 1,2 / 3,4 |
| 1 | 8289654 / 8289655 | 8289642 / 8289643 | 5,6 / 7,8 |
| 2 | 8289656 / 8289657 | 8289644 / 8289645 | 9,10 / 11,12 |
| 4 | 8289658 / 8289659 | 8289646 / 8289647 | 13,14 / 15,16 |
| 8 | 8289660 / 8289661 | 8289648 / 8289649 | 17,18 / 19,20 |
| 12 | 8289662 / 8289663 | 8289650 / 8289651 | 21,22 / 23,24 |

Each accession has the prefix GSM. Every one of the24 processed column names matches its official GSM title after replacing spaces with underscores. `sample_map.json` and`.csv` retain all28 full accessions, BioSample/SRA experiment links, platform and processed-column coordinates.

The additional GSE268324 samples are CNOT3-IP-1=GSM8289671, CNOT3-IP-2=GSM8289672, Input-1=GSM8289673 and Input-2=GSM8289674. The paper calls selective-ribosome-profiling replicates biological; this does not establish independent SLAM preparations, exact input/IP culture pairing, cross-condition pairs or time-course culture linkage. The assays share a cell-line/study context, not authenticated common biological preparations.

## What the processed endpoints actually contain

| Endpoint or requirement | Established | Unresolved limit and eligibility consequence |
| --- | --- | --- |
| Processed trajectories | Original workbook columnsC:Z are explicitly documented as transcript T>C conversion ratios,24 time/condition/replicate-labelled columns; Gene_Name and Gencode Gene accession accompany them. | Fraction-versus-percent convention, numerator/denominator detail, background correction and spike-in normalization are not specified by these headers/legend. No coverage/count/error columns occur in the inspected header schema. |
| Source processing | GEO states GRCh38 alignment and slamdunk0.3.3; other assay processing mentions HISAT2.1.0 for ribosome reads. | These declarations do not specify a half-life estimator, replicate aggregation or its quality/censoring rules. |
| Paired reported half-lives | Correct TableS2 has `HL_WT` and`HL_CNOT3KO`, described as half-lives from Zhu2024, with one field for each condition per GeneID. | They are not labelled replicate-specific, and the metadata do not say arithmetic mean, pooled fit or another aggregation. Do not call them averages. Missingness/valid-positive denominator has not been examined. |
| Relative half-life response | `HL_FC_LOG2` is explicitly log2(HL_CNOT3KO/HL_WT). | A frozen audit can verify this identity and compare ranks after outcome authorization. Units cancel only if the intended common scale is respected. |
| Absolute rate response | No calibrated rate field occurs. Reciprocal half-lives can be calculated algebraically under a specified convention. | Stored HL time unit remains unauthenticated. Chase times being in hours does not prove these stored fields are in hours. No h^-1 claim is warranted. Use a clearly stipulated common deposited-time unit only, or defer an absolute-rate endpoint. |
| Measured pathway association | `LOG2_FC_IP`, `STAT_IP`, `ADJP_FC_IP` and`SE_IP` describe CNOT3-IP versus regular ribosome profiling using DESeq2. | `SE_IP` is association fold-change uncertainty, **not half-life SE**. It cannot propagate uncertainty through reciprocal HL or provide a decay confidence region. The method/caption discrepancy about DESeq2 versus Fisher must remain visible in reproduction. |
| Half-life uncertainty/quality | No SE_HL, replicate-HL, confidence interval, fit R2, cap/censoring or explicit half-life-fit-quality field occurs in the TableS2 header/description. | These are absent from the authenticated schema, not merely unread numeric values. Cited Herzog-method rules cannot be imported as authenticated Zhu settings. |
| Growth dilution | Primary-prose review reports an original supplementary growth effect under perturbation. | Per-condition growth correction in the deposited HL estimator is not authenticated. Thus reported decay/stability change is not identified as molecular degradation change alone. |
| Perturbation specificity | Source labels/control and CNOT3-targeting guides are traceable. | Residual activity, compensation, cell-state change and other decay pathways are not ruled out by a pair of reported HL numbers. A causal pathway contribution cannot be equated to their difference. |

The physical worksheet dimension is not an outcome-valid sample size. TableS2 contains formatted rangeA1:X22849; original trajectory sheet rangeA1:AH19799, but only the named26 header fields are documented. No claim is made about unlabelled/uninspected cells elsewhere in those sheets.

## Identifiers and the reference-feature gate

Identifier-only inspection, explicitly authorized by root, found:

- Original trajectory workbook:19798 unique GeneIDs, no duplicates;19780 are versioned ENSG identifiers and18 have other forms.
- Reanalysis TableS2:22848 unique unversioned ENSG identifiers, no duplicates. Inventory `reanalysis_gene_ids_only.json` SHA256`85fb60552f3c12477d28dbeeafbd171d52828ecb66e30591cefa60f77f482c52` was supplied to the reference-only lane. It contains row/cell/type/ID only.
- Raw exact identifier overlap is zero. A separately disclosed diagnostic strips only a terminal dot+integer from an exactly matching ENSG identifier, with zero collisions. It gives19563 common gene cores,217 original-only cores and3285 reanalysis-only IDs;18 other original forms remain unmapped. This is a name-level join, not proof of common source annotation, transcript sequence or expressed isoform.

GEO declares GRCh38 but no GENCODE/Ensembl release, GTF filename or exact annotation download. The original legend calls IDs Gencode Gene accessions. The reanalysis methods, as authenticated by the prose lane/root, describe Ensembl GRCh38 longest valid CDS without a release. TableS2 has no transcript identifier, sequence or full codon-composition columns. `GENE_LEN` is described as longest possible transcript length, so it must not silently be relabelled measured CDS length.

**A new pinned representative-CDS feature set is defensible for a conditional gene-level predictive question.** Select it by a deterministic annotation-only rule, authenticate release/checksums and gene-ID mapping, preserve all failures and expose the selected transcript/sequence/feature provenance before viewing target values. It is a newly defined reference feature set, not an exact reproduction of the unspecified source release or the codon composition of the expressed molecules. Association incremental to those reference features does not prove incrementality beyond true expressed-isoform composition or all UTR/coverage confounding. That stronger transcript-specific inference gate fails. No sequence/reference acquisition was performed in this audit.

## Smallest defensible next analysis and reproduction prerequisite

Freeze an all-eligible-GeneID point-measurement audit before opening HL or IP values. Predeclare missing/finite/positive-HL rules, identifier joins, annotation mapping, tie handling, quality fields, exclusions, the primary response and all challenges. Do not select the most favorable genes, quantiles or endpoint after inspection. Preserve negative reciprocal-HL contrasts rather than clamping them.

The safest unit-independent endpoint is the explicitly defined log2(KO/WT) half-life contrast. A descriptive reciprocal contrast `ln2/HL_WT-ln2/HL_CNOT3KO` and its ranks can also be examined **under an explicit assumption of one common positive deposited-time scale**, labelled accordingly. Global positive rescaling changes its magnitude but preserves sign/rank. Condition- or gene-specific scaling, undocumented censoring or different growth correction do not enjoy that guarantee. No hour-calibrated rate, causal pathway fraction, experimental confidence interval or biological noise decomposition follows. Baseline HL appears algebraically in the response, so mathematical coupling and measurement error require explicit challenge; a baseline-adjusted predictive coefficient is not a causal correction.

The source-backed reanalysis Figure2 reproduction definition, supplied by the separate prose lane, is: genes with at least10 aligned reads in each of2input and2IP libraries, with reported half-life availability(N=5481); rank DESeq2`STAT_IP` descending; top1000 versus remaining4481, with no adjusted-P threshold for that contrast. TableS2 does not expose the four raw read counts or an explicit coverage-pass flag, so reproducing that threshold from this table alone is not yet authenticated. Complete-case reconstruction must first be checked against the published denominator. Exact missingness/tie rules are not specified. Any top-PTMD comparison has an additional unresolved question about inherited versus reranked top1000 membership. These are reproduction obligations, not grounds to substitute a favorable cohort.

Only after that control reproduction and a frozen design should a gene-held-out comparison test additional association information beyond the **pinned reference** full codon/amino-acid composition, length/GC, baseline and available IP-quality proxies. Gene holdouts assess cross-gene prediction in this context, not independent experiment validation. Missing half-life error/coverage data remain a limitation, not a field that can be imputed from SE_IP. A trajectory refit would be a different prespecified descriptive measurement model, not a silent source-fit replacement.

## Exposure and stopping record

The exposure ledger distinguishes metadata/documentation from outcomes. Published cohort sizes, top1000 selection,17 significant IP genes, overlap279 and published correlations about0.18/PTMD and0.06/IP were communicated by the primary-prose lane. They are prior published aggregate exposure, not calculations from the protected gene-level table, and the proposed test is not blinded. Searching the already-extracted supplementary prose also exposed growth/replication captions; these were not transferred to the SLAM replicate labels.

No data-bearing row's HL, association, conversion ratio or quality value was extracted. No actual valid-positive or complete-case denominator was computed. No calibration or validation cohort was selected. No raw pipeline, new hypothesis, Git action, publication or shared-file edit occurred. This audit stops at this eligibility decision. Specific source URLs/hashes, all28 samples and the precise field-access ledger are retained in this owned directory.
