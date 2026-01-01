# cPRIME B0: marker-context and genotype identification gate

**Decision:** qualified metadata GO for an observed-local-sequence marker-context audit, but **not yet freeze-ready for an isolated intended-edit effect, a complete factorial interaction, or full realised-haplotype claims**. Matched marker contexts genuinely exist. However, the source already demonstrates their non-neutrality and compares deliberately confounding markers. A new analysis must target a precise invariance or identification question rather than repeat that result. One arrayed genotype label has a consequential unresolved sequence mismatch.

No variant phenotype, count, RPM, score, outcome missingness, notebook output or raw read was inspected. The only outcomes exposed are already published statements/aggregates in primary prose and an incidental search snippet, logged below. No modelling or phenotype-based selection occurred.

## Charter

Goal: decide whether observed intended-edit/marker states identify the proposed effect. Scope/object/subject: retrospective endogenous RNA measurement / cPRIME NF1 local sequence states / genotype, context, amplicon and preparation matching before outcomes. Type: descriptive source-contract audit. Questions: Q1 what exact object is counted; Q2 which matched states exist; Q3 what replication/context/transcript information is retained; Q4 what distinct, executable inference remains. Sources rank the originating paper and supplements first, pinned author code/metadata second, then official accession and bibliographic records. The 20-page expansion plan's pages 1-11 and full latest user request were read.

## Canonical sources

| Rank | Source | Authentication and use |
| --- | --- | --- |
| 1 | El Kassem et al., cPRIME, Nature Communications 17:10872 (2026), DOI [10.1038/s41467-026-78344-9](https://doi.org/10.1038/s41467-026-78344-9), PMID 42888231, PMC13675207 | Europe PMC full-text XML, DOI/title matched; Methods, marker section and selected replication prose inspected. Published 15 October 2026, accepted 22 September 2026. Publisher indexed page still describes early accepted publication subject to editing, so do not silently assume immutable final wording. |
| 1 | Official supplementary inventory and Data 1-9, Source Data 1/2 | Retrieved through the public PMC cloud inventory. Explicit file mappings, schemas and whitelisted metadata only. Data 6 is arrayed validation, Data 7 is exon-27 design, Data 8 is multi-exon design/measurements/primers. |
| 1 | [Boettcher-Lab/cPRIME](https://github.com/Boettcher-Lab/cPRIME/tree/30125415304664c33370cd7801f820a345834e40) | Commit/tree 30125415304664c33370cd7801f820a345834e40; downloaded only relevant metadata/search sequences/counting code and source-only analysis cells. No author code executed. |
| 1 | [Zenodo 21991157](https://doi.org/10.5281/zenodo.21991157) | Record metadata authenticates v1.0.0 archive, 45,104,767 bytes, deposited MD5 c2a125f3d630703032100a8bd548b948. Large archive not downloaded because the pinned repository supplies this gate. Do not claim archive-byte equivalence verified. |
| 1 | [PRJNA1257064](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1257064) | Official ENA file-report metadata: 330 run/library rows, no reads or read counts. This is an archive inventory, not 330 independent preparations. |
| 2 | Crossref and Scholar Sidekick correction check | DOI record matched; no indexed retraction, correction or concern returned on 25 October 2026. This is a checked-record status, not a guarantee that no future notice exists. |

Detailed URLs and byte hashes are in the acquisition receipts and final manifest. Nature browser fetch failed through an inaccessible identity redirect; PMC browser returned a challenge. Neither was bypassed. Normal public Europe PMC XML and independently documented public PMC cloud objects succeeded. The Europe PMC supplementary-files API timed out and the Springer CDN hostname failed DNS; both failures are retained, not reported as permission denials.

## Findings ledger

| ID | Source-supported finding | Evidence | Answers |
| --- | --- | --- | --- |
| F1 | Exon-27 design has 98 library rows and 16 additional author validation rows. The published validation sheet has 13 rows; these inventories are not silently substituted. | Author Overview workbook; published Data 7, whitelisted sequence/ID cells | Q2 |
| F2 | Multi-exon design has 341 IDs and 178 intended-edit labels. Both DNA and RNA search sets have exactly six duplicate-motif pairs: five legacy `_old` aliases and one reversed intended-edit/marker pair. Thus they contain 335 distinct motifs. | All 12 pinned multi-exon search-sequence text files; source-only analysis deduplication lines 30-38 | Q1,Q2 |
| F3 | Data 8 contains 335 shared combined IDs in each of KBM7 and K562, three metadata rows per ID/context: 2,010 rows. Its 177 intended labels include 111 with multiple retained marker states and 66 with one. | Only Data 8 columns 1-7, 27 and 33; matched author design IDs | Q2,Q3 |
| F4 | The reference-annotated synonymous subset has 63 retained intended labels when `splice_region_variant&synonymous_variant` is included: 33 have one marker state, 24 two, six three. Thus 30 have multiple states. | Annotation strings and canonical design IDs, no phenotype field | Q2 |
| F5 | The R1207R arrayed panel lists intended-edit-only, marker-only and alternative combined-marker states. The row labelled `marker_GtoT` has the **same DNA and RNA motif as the original G>A combined state**, in both author text files/workbook and published Data 7. | `published_validation_sequence_metadata.json`; exact motif comparisons in contract receipt | Q1,Q2 |
| F6 | What is counted is a short exact search motif containing edit plus marker and nearby WT flanks, not an independently called complete realised genomic haplotype. | Methods Sequencing and data analysis; count_ss.py; Count_edits.sh | Q1 |
| F7 | cDNA detection requires preserved search sequence, flanking primer sites and an approximately canonical gel-selected product. The detected-isoform signal is not automatically total RNA. | Primary explanation accompanying Fig. 2; cDNA amplification Methods | Q1,Q3 |
| F8 | The source labels n=3 biological replicates. Initial screens are explicitly transduced in triplicate; arrayed Methods also mention frozen triplicate aliquots but do not fully specify independence of transduction/culture batches for every condition. | Main Figs. 1-3, supplementary Fig. 1 replication prose; Methods | Q3 |
| F9 | Data 8 omits an explicit replicate field even though analysis pivots on cell line/exon/replicate. ENA has all 72 unique multi-exon names DNA/RNA × six exons × two cell lines × R1/R2/R3. | Header schema; analysis source lines 59-77; ENA metadata | Q3 |
| F10 | The source already compares markers, removes predicted-confounding markers before intended-edit averaging, and reports arrayed marker effects. | Primary marker section; SNV class assignment Methods; source-only analysis | Q4 |

F2 resolves the 341-versus-335 inventory difference without reading coverage or outcomes. The reversed-role duplicate also explains why 178 intended labels become 177 retained labels. Do not treat two guide/design labels for that same motif as independent genotype measurements. F4 counts reference annotations, not measured protein preservation or an empirical phenotype-valid denominator.

## Exact object and unit contract

| Link | Authenticated object | What remains unresolved |
| --- | --- | --- |
| Intended edit → marker | Data 7/8 have separate genomic IDs, HGVS annotations, combined IDs and ENST00000358273 annotation. Exon-27 design uses hg19/GRCh37; multi-exon uses GRCh38. | Do not join raw coordinates across builds or legacy labels without the recorded mapping. |
| Marker → measured genotype | Local WT-versus-edited search words authenticate the named substitutions within their flanks. DNA/RNA motifs match for the examined arrayed R1207R states. | Additional changes outside the word, distant phase, alternative alleles and full realised haplotypes are not excluded by a substring match. Counting can match motifs, rather than mutually exclusive haplotype classes. |
| Genotype → cell context | All 335 canonical multi-exon motif IDs occur as metadata in both cell contexts; each has three rows. | Cell line bundles ploidy, clone/background and other context. Initial exon-27 harvests differ: KBM7 day 21 versus K562 day 45. Multi-exon harvest is day 45. A cross-cell contrast is not a pure single-factor context intervention. |
| Cell sample → DNA/RNA | Methods derive gDNA and RNA from source cell preparations; file names encode modality/context/replicate. | Data 8 does not preserve explicit per-row replicate keys or sample-accession joins. Three row positions are not automatically authenticated biological pairings. Seventy-two amplicon libraries are not 72 independent cultures; exons share source material. |
| Amplicon → isoform | Exon-27 gDNA primers flank introns 26/27; cDNA primers lie in exons 26/28. Data 8 supplies per-exon primer metadata. | Reference transcript annotation is not detected-isoform assignment. No complete per-genotype junction/isoform census is authenticated by these processed files. Excluded/size-shifted transcripts can be invisible. |
| Counts → endpoint | Paper endpoint is cDNA RPM/gDNA RPM for a local variant-marker word. Multi-exon analysis uses `(count+1)/(total_reads+1)*10^6` then ratio, and deduplicates identical words before joining design metadata. | Data 8 exposes RPM/score schema, without explicit raw count/denominator or replicate columns. Freeze source pseudocount/denominator rules before an uncertainty model. Published descriptions of counts do not prove every processed sheet contains raw counts. |

The archived data provide useful local sequence-conditioned responses. They do not establish a pure intended edit merely because a guide was designed for it. Likewise, an unchanged reference protein annotation is not measured preservation of all processed products.

## R1207R state inventory and blocker

Reference annotation: intended c.3621A>G; original marker c.3624G>A. Sequence comparison, rather than labels alone, gives:

| State | Deposited ID/label | Local substitutions actually encoded by the supplied motif |
| --- | --- | --- |
| WT | A WT search word is supplied as design metadata | No mutation. An explicit matched WT cDNA/gDNA response is not authenticated in Data 6; it cannot simply be set to one. |
| Edit only | `Synonymous_edit_29560144_wo_marker` | c.3621A>G |
| Marker only | `Synonymous_edit_29560144_wo_intended_edit` | c.3624G>A |
| Combined A marker | `Synonymous_edit_29560144` | c.3621A>G plus c.3624G>A |
| Combined C marker | `Synonymous_edit_29560144_marker_GtoC` | c.3621A>G plus c.3624G>C |
| Labelled combined T marker | `Synonymous_edit_29560144_marker_GtoT` | Supplied motif instead encodes c.3621A>G plus c.3624G>A, identical to combined A marker |

The last row is a metadata/count-reference discrepancy, **not proof that the physical edited cells or the paper's biological observation were wrong**. No outcome can repair the sequence identity by itself. Do not silently replace its A with T or treat the labels as three distinct realised contexts. The relevant source record to close it is the actual counted motif/sequence-verified state and sample mapping for that condition, not its phenotype.

## Scope filter

Keep local-sequence marker-context contrasts separate from isolated intended-edit effects, full haplotype effects, total RNA, protein function and clinical interpretation. No new editing designs or experiments enter. Published marker non-neutrality is prior, not this project's discovery. All primary marker examples are already exposed through the paper and cannot be described as untouched biological validation.

## Object filter

Canonicalize duplicate DNA/RNA words before counting contexts. Preserve distinct assemblies and incomplete annotation. Treat source-labelled biological replication as such, while retaining undocumented independence/pairing rather than inferring it from R1-R3 suffixes. Avoid using the six metadata aliases as outcome-dropout evidence. Do not infer true zeros or observed coverage from any metadata inventory.

## Conclusions

1. **Combined-state and marker-choice effects are potentially measurable** for canonical local states in a fixed cell/amplicon, once replicate keys, source normalization and an uncertainty contract are frozen (F2-F9). At least 30 reference-synonymous intended labels have multiple marker states in both cell contexts. This is a finite metadata eligibility count, not a prevalence estimate or result.
2. **Intended-edit effects are generally not identified by those combined states.** On a declared response scale let Y_em=mu+alpha_e+beta_m+eta_em, with baseline conventions at e=0,m=0. Without edit-only/WT or restrictions on both beta and eta, any alpha can be offset by eta while matching all observed combined states. Marker-only bounds do not restrict alpha when interaction is unrestricted. Even a measured Y_e1-Y_01 equals the intended effect plus interaction, not the isolated intended effect. If both |beta|≤B and |eta|≤I are justified on that scale, the elementary interval around Y_e1-Y_00 has half-width B+I; these bounds are assumptions, not estimates available from this gate.
3. **A four-state interaction contrast requires a measured, comparable WT response too.** The one marker-only panel is promising but the WT endpoint and the T-marker identity are unresolved (F5). RPM-based cPRIME scores are not automatically WT-normalized effect sizes.
4. **Broader context portability is a conditional reserve.** The 335 exact local motif IDs overlap across cell lines, but transcript visibility, duration/background and missing replicate keys limit attribution (F7-F9). No full-haplotype or transcript-assignment claim follows.
5. **Novelty is narrow.** The original source already establishes marker effects and removes selected confounding markers. A useful increment would be an all-declared-context invariance/identification audit with a fixed practical decision rule, preserved uncertainty and complete intended-edit holdout. An additional marker-difference plot or the selected R1207R example alone repeats the source (F10).

## Next steps / Risks

**Smallest defensible prospective contract:** first freeze a metadata-only join of canonical combined word → intended edit → marker → cell/amplicon → actual DNA/RNA sample and biological replicate. Resolve the R1207R label discrepancy before including that state, and authenticate whether a matched WT response exists before any factorial claim. No raw-read reprocessing is part of the present gate.

If those records are recovered, a low-cost empirical audit can retain all 30 reference-synonymous multi-context intended labels, with all their declared marker states and source-supported replicates, and ask whether a fixed observed-local-state conclusion survives marker choice in each cell context. Assign any development/validation division by whole intended edit **before outcomes**, preserving shared-marker/genotype dependencies. Fix endpoint, practical tolerance, all-state missingness/zero rules, pseudocounts, reference and covariance procedure first. Do not remove predicted bad markers and then claim invariance over the original marker set. Such a result concerns the **observed finite marker set**; it does not isolate alpha or establish all-marker invariance.

No empirical opening is authorized by this report. If the needed sample/replicate join or comparable WT state cannot be authenticated, retain only a bounded combined-state descriptive question and the explicit nonidentification statement. The selected arrayed panel remains source-exposed validation, not new external validation. Untouched whole intended-edit groups are a computational exposure boundary within the same source, not independent experimental replication.

## Exposure and acquisition record

- Inspected all requested plan pages 1-11 and full current user request; no plan claim was accepted without a source check.
- Read primary Methods and marker-context prose, including published named examples, their qualitative phenotype conclusions, reported nonsense averages 0.19/0.30 and Pangolin example scores. An incidental third-party search snippet also displayed the same nonsense averages; it was not used as evidence. These published exposures are not newly analyzed data and cannot be called untouched.
- Read first-row textual headers/sheet names/dimensions only from mixed source workbooks; numeric first-row cells were masked. Sheet dimensions are not phenotype-valid counts.
- Read only explicitly whitelisted sequence/ID/annotation/context columns from Data 7/8, the Data 6 ID/name columns, and design-only author tables. Did not inspect count, RPM, cPRIME, class, prediction-score cells or their missingness. `select_metadata.py` records the whitelist and does not return phenotype columns.
- Downloaded the source notebook but inspected **code cells only**, never outputs or embedded figures. Code was not executed. Its source normalization constants and paths were read as methods, not outcomes.
- Retrieved ENA sample/library metadata only. No raw reads, counts, fitting, author contact, credentials, external writes, Git or publication actions.
- Primary XML/PDF/workbooks and the output-bearing notebook remain local evidence cache. They are not a whitelist for redistribution in the computational archive. Promote reports, receipts, reader and necessary selected metadata only after normal archive review.
