# PERSIST C4 measurement eligibility, cycle 3

## Decision

A bounded **source-aware cellular-decay prediction comparison** is feasible in one exact NanoLuc protein and transcript-context stratum, conditional on the response-scale check below. It is not an independent-study validation, a complete synonymous landscape, a measured structural-barrier test, or a library-wide measured-protein comparison. No cellular response values were inspected and no prediction was evaluated. The proposed freeze requires the root agent's review before outcomes are opened.

The strongest internally consistent candidate has 37 deposited processing-eligible rows, each encoding the same 221 amino acid NanoLuc product, sharing transcript sequence outside the CDS after replacing the unique experimental barcode with a placeholder. There are 22 prefix-310 and 15 prefix-930 designs. Their source labels and ancestry must remain explicit. This is one study and one protein, with two design-provenance strata, not 37 independent biological families.

## Charter and ranked sources

This descriptive audit asks which already-measured constructs support a permitted, identifiable endpoint and validation unit. Scope: retrospective nonpathogenic reporter measurement. Object: the full-length PERSIST reporter library. Subjects: deposited construct sequences, identifiers, context, sample labels and processing definitions. The knowledge depth was working: an earlier audit had established endpoint types but not exact eligible denominators.

Rank 1: [Leppek et al. 2022, DOI 10.1038/s41467-022-28776-w](https://doi.org/10.1038/s41467-022-28776-w), especially Methods and replication statements. Public JATS was recovered from [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8940940/fullTextXML). Rank 1/3: [author repository](https://github.com/barnalab/persist), exact blobs in tree `8922dc115acf0c3c803b15fbc97e67edce805c36`. Rank 3: reproducible inspection of metadata and identifier projections. Earlier local audit/manifest provide context, not independent corroboration. No secondary-source conclusion is needed.

The PMC webpage returned a browser/reCAPTCHA page. No challenge was attempted. The ordinary public Europe PMC XML request succeeded, as did GitHub blob requests. Receipts, content hashes and pinned blob identities are in `source_manifest.json` and `retrieval_receipts.json`. There is no source-unavailable claim.

## Exact inclusion accounting

The starting table has 233 rows. A conservative, disjoint metadata filter gives:

| Step | Rows |
| --- | ---: |
| Exclude MEV or CoV-labelled construct groups | 29 |
| Exclude remaining explicitly viral-derived UTR rows | 47 |
| Exclude remaining `others` or unresolved UTR-origin rows | 38 |
| Retain both UTRs labelled cellular, with nonexcluded reporter groups | 119 |
| Retained rows with deposited `longamp=TRUE` | 105 |
| Retained rows listed in deposited solution-output identifiers | 105 |

This conservative filter leaves some potentially benign synthetic UTRs unused because their ancestry is not closed. It does not reclassify a design as pathogenic because its informal name contains terms such as “OpenVaccine” or “viral letter frequency”. Inclusion here is restricted to the actual NanoLuc/eGFP reporter, cellular UTR labels and verified coding sequence. Vaccine/MEV and CoV construct groups remain excluded and were not translated or evaluated.

Across the complete deposited table, `longamp=TRUE` has 192 rows. The headerless solution output also contains 192 unique SequenceIDs, and the two sets are exactly equal. The MATLAB source comment says the final table has 191 entries; that documentation count disagrees with the actual deposited file. The older “203 minus 12 equals 191” account must therefore not be used as a verified deposited denominator. No outcome-based re-filtering was attempted.

The exact set equality is also a provenance warning: `longamp` has not been established as a sequence-only coverage criterion independent of the study's outcome-based processing. The original solution analysis filters anomalous time courses. Use this flag as an **author-deposited processing eligibility flag**, disclose potential selection dependence, and do not claim that the retained rows are an unbiased sample of all designs. Presence of a construct identifier does not establish a finite/nonmissing cellular coefficient, measurement quality, or calibrated uncertainty; those checks remain gated with outcome access.

## Exact construct and protein identity

Every one of the 119 retained rows has a nonempty `CDS.sequence.longest`, uniquely located in `RNA.sequence`; the reported full RNA length agrees with sequence length. Each CDS has an initial ATG, an in-frame terminal stop, no internal stop, and agrees under two translation routes: Biopython and an independently implemented standard-code lookup. This establishes the identity of the deposited intended construct, not resequencing confirmation of physical molecules.

Twenty-four rows have an empty `CDS.sequence` and an 850 nt `X5UTR.sequence`, despite nominal hHBB labels. The longest-CDS field and full RNA sequence recover their coding identity. Treating the empty CDS field as an empty protein gives a false family. The audit script corrects that interpretation without editing the deposited data.

The retained set contains 21 exact translated-protein identities. The main repeated families are:

| Exact protein and context | All rows | Distinct CDSs | `longamp=TRUE` | Solution-ID overlap | Joint flagged/listed |
| --- | ---: | ---: | ---: | ---: | ---: |
| 221 aa NanoLuc, context A | 21 | 17 | 21 | 21 | 21 |
| Same 221 aa NanoLuc, context B | 43 | 43 | 37 | 37 | 37 |
| 284 aa eGFP-associated product, common context | 9 | 9 | 8 | 8 | 8 |
| 264 aa eGFP-associated product, common context | 2 | 2 | 2 | 2 | 2 |

Other eGFP-labelled rows encode three singleton protein identities. Thirteen additional NanoLuc-design rows each encode a different exact 221 aa product; they must not be pooled into the synonymous main family. There are also 23 rows with the same 206 aa starting NanoLuc product in 23 distinct contexts, and five rows with a 212 aa product in five contexts. “Nluc”, “NlucP”, “eGFP”, and a shared UTR name do not establish exact protein preservation.

The main 221 aa protein SHA256 is `dd0e8ca4f9a41c1fdd610ac0e48d953ff55be9483e2067db5e713804931fc148`. Its barcode-normalized context hashes are A `449b0eb3d167c8cdb1f8e6f36bf6d3961b46a8bd6da9310c358ac35a2ba6e6f8` and B `e682ac5099bf1b6f4beb22f1f9817c89022c0bd31b279975a872546eddf726e9`. Both transcripts are 928 nt, but the barcode/constant cassette occupies a different location relative to the downstream UTR sequence. Shared `hHBB` names therefore do not justify pooling A and B as one identical context. Exact identifiers are in `eligibility_results.json`.

Context A contains 16 prefix-310/Eterna-labelled rows, four RiboTree-labelled rows and one MCsampler-labelled row. Four excess barcode-associated rows reduce its 21 rows to 17 distinct CDSs; barcode variants cannot be split across training and validation. Context B's 37 flagged rows comprise 22 prefix-310/Eterna-labelled and 15 prefix-930/RiboTree-labelled rows. The minimum cross-stratum CDS Hamming distance is 118. These labels record deposited provenance; a crowd submission may itself use an algorithm, and neither stratum is an independent study or necessarily an independent design process.

No pair of full RNA sequences in any repeated fixed-protein/context stratum differs by one nucleotide after retaining its actual experimental barcode. Ignoring barcodes, all 43 context-B CDSs have five one-nucleotide edges; the 37 flagged CDSs have four. The latter consists of two three-vertex paths and 31 isolated vertices. Context A and both repeated eGFP strata have no one-nucleotide CDS edges. Thus this source does not furnish a complete measured local mutation landscape, even though it contains some isolated CDS-neighbour comparisons. Replacing a barcode with a placeholder is a grouping convention, not a claim that the measured RNAs were otherwise identical genotypes.

## Endpoint and replication interpretation

| Endpoint | Actual meaning and available representation | Eligibility limit |
| --- | --- | --- |
| Cellular persistence | Time-course-derived decay of RT-amplifiable reporter RNA in HEK293T; stored coefficient, coefficient error and half-life fields | A fitted RNA persistence proxy, not cellular fitness or protein amount; no individual coefficient values inspected |
| Solution stability | Condition-specific RT-amplifiable RNA decay in a separate solution assay; identifiers and fitted-output schema available | Different chemical/capping/poly(A) context from cellular measurement; not interchangeable with cellular persistence |
| Ribosome load | Weighted normalized RNA abundance across polysome fractions | A translation-associated proxy, not directly measured protein production |
| Predicted expression | Kinetic calculation combining RNA decay, ribosome load, molecular mass and protein-stability assumptions | Derived from other endpoints; cannot serve as independent measured protein validation |
| Measured protein output | Selected NanoLuc/Firefly luciferase follow-up panels | Not a label for all 233 rows; source Methods state signal peptides were removed for these assays, so match the actually assayed construct before joining |

The sample sheet has three Long-labelled records at each of 1, 7, 12, 24 and 106 hours, and matched Short amplicon records. These are not six independent experiments. The author R script selects only Long samples at 1, 7 and 12 hours for cellular fitting: Sample_IDs 13,14,15,18,19,20,27,28,29. It estimates unmoderated regression coefficient errors. The sample sheet does not identify the biological origin of the repeated time series or provide culture/extraction pairing labels. Do not relabel these errors as independent biological-repeat variability or infer a biological noise ceiling.

The polysome sheet has a spike-in and a no-spike-in fraction series; processing uses the spike-in series. Those normalization variants/fractions are not independent biological replicates. The published solution half-life caption explicitly describes one biological replicate with bootstrap uncertainty. The 100-resampling setting in the MATLAB source resamples timepoints, not independent biological cultures. Selected luciferase assays have their own replication statements, including technical triplicates in the formulated-RNA panel. None supplies missing library-wide biological replicate identity.

The main-library methods distinguish the cellular capped/polyadenylated pool from the solution pool. Chemical modification is not a per-row factor in the processed construct or sample table. Restrict a proposed analysis to the documented main cellular pool; do not pool selected pseudouridine/formulated follow-ups or infer modification status from the label `NlucP`.

**Response-scale gate:** the R source fits `log2(counts)` versus time and sets the stored decay slope to minus the fitted coefficient. It then calculates `log(2)/degrate_ics`. A negative log2-abundance slope has half-life `1/slope`, whereas a natural-log slope has half-life `ln(2)/slope`. This is an algebraic source-code inconsistency, not an inspected outcome discrepancy. Use a clearly specified slope endpoint and verify how the released column was generated before interpreting any half-life or physical decay-rate units. No creator code was executed.

## Proposed one-endpoint freeze, not executed

The narrow question is whether sequence features improve prediction of **cellular RNA-decay slope across deposited design-provenance strata**, for context-B main-protein constructs. It does not test a newly discovered RNA mechanism or mutation-barrier theorem.

1. Freeze the 37 explicit context-B identifiers, exact protein/context hashes, and the deposited selection flag. Before opening their values, resolve the slope-column scale and document treatment of nonfinite or missing responses without replacing them by zero. Keep all 37 in the intended denominator and report any later measurement exclusions. Do not select using predicted expression or solution outcome values.
2. Use the 22 prefix-310 rows as development and hold the 15 prefix-930 rows for one transfer evaluation, or reverse this direction only if fixed before outcome access. Do not choose the direction from performance. Report these as design-source strata within one reporter/study. Keep identical/near-neighbour CDS groups together in any internal validation.
3. Primary endpoint: the source-defined negative log2-abundance slope, or a documented proportional natural-log conversion fixed before fitting. Do not optimize transformations. Primary comparison: absolute prediction error on all eligible held-out rows, reported per row and as a mean without an independent-gene confidence interval.
4. Serious baselines: training-set mean; a low-dimensional composition model with GC and prespecified nucleotide/dinucleotide summaries; and codon-composition ridge regression with tuning confined to grouped development data. Length is constant in this stratum and cannot explain within-stratum differences, but should remain in cross-context bookkeeping. A source label alone is constant in a wholly held-out source and cannot estimate its new intercept; use a training-only fallback rather than calibrating on validation labels.
5. A structural feature model would need a separately frozen, low-complexity specification and the same split. Released stability/structure predictors are comparison candidates, not new independent response data. No feature set, hyperparameter or threshold may be selected using the 15 held-out responses. Preserve negative transfer if it occurs.
6. Stop if the scale cannot be established, too many held-out outcomes are unavailable, the design matrices are not estimable under the frozen model, or the only apparent improvement requires validation recentering. Do not rescue failure by pooling context A, altering protein identity or switching to model-derived expression.

The 8 flagged 284 aa eGFP rows offer a possible later protein-family transfer check, but they change protein and transcript context jointly and provide only one small second family. The two 264 aa rows are not an adequate independent validation family. No claim of clean protein-only transfer or broad biological generalization is warranted. Context A offers another source/context sensitivity check, not another gene or study. The 23 starting-NanoLuc UTR contexts test a different question with a fixed CDS. Keep it separate.

## Findings ledger and filters

| Finding | Evidence | Interpretation |
| --- | --- | --- |
| 119 conservative rows; 105 flagged/listed | Pinned `constructs.tsv`, projected metadata, solution identifier projection; `eligibility_results.json` | Exact metadata denominator, not inspected outcome completeness |
| Protein and context labels hide distinct constructs | Unique CDS-in-full-RNA check, two translations, barcode-normalized context hashes | Enforce exact identity, not group names |
| 192 deposited solution rows versus 191 comment | Identifier-only row count and MATLAB source comment | Documentation discrepancy remains explicit |
| Source-aware 22/15 split exists | Fixed context-B identifiers and prefixes, name metadata; minimum cross-stratum distance118 | Feasible within-source transfer comparison, no independent-study unit |
| Response scale requires reconciliation | `analysis/persist_p.R`, `log2`, `ICS_USED`, `degrate_ics`, `halflife_ics` assignments | Freeze a verified slope endpoint before outcomes |
| Replicate independence not established | Sample sheet lacks biological-origin column; Methods and caption distinctions | No biological-noise calibration or unsupported independent-repeat claim |

Scope filter: excluded viral-derived UTRs, vaccine/MEV and CoV groups, uncertain UTR-origin rows, intervention design and outcome evaluation. Object filter: revised nominal protein families into exact sequence identities; separated contexts, barcodes, endpoints and source strata; rejected a complete-landscape or independent-study interpretation. The remaining claims stay within retrospective reporter measurement eligibility.

Research progress: canon and charter recorded; descriptive/source instructions read; ledger and both filters completed; conclusion supported by metadata checks, with response scale, selection provenance, output missingness and biological-repeat identity left as explicit gates.

## Exposure and execution record

The previous audit already encountered published aggregate findings and selected named reporter highlights. This cycle read Methods, metadata and caption sentences concerning replication. A replication-caption extraction also displayed published significance notation and the text `n=3`; these are aggregate source-text exposures, not a new evaluation. An overly broad initial local search returned snippets from unrelated earlier literature records; those were not pursued or used and no excluded analysis branch was reopened.

**Accidental solution-row exposure:** the solution output has no header. The initial retrieval printed its first line under the mistaken label “header”, exposing the two fitted numerical fields for SequenceID `130006D1`, named `humanHBB_ModernaRef_egfp_lowstructure_HBB_F30Pepper`. Their values are intentionally not repeated here. This specific solution-stability row is not untouched. The incident was immediately reported to the root; no fit, comparison, threshold or decision used those values. Subsequent extraction retained only the first three identifier/name fields. No cellular coefficient/half-life, ribosome-load, predicted-expression or protein-outcome cell was displayed, numerically inspected or evaluated.

The combined deposited construct table and solution table were downloaded as bytes for preservation and field projection; this is not a claim that their raw files contain no outcomes. Projected metadata files are the only inputs read by `eligibility_metadata_audit.py`. The whitelist excludes all response and precomputed predictor columns. The audit executes sequence identity, translation, context, identifier-overlap and sample-count checks only. Author R/MATLAB code was read but never executed. Reviewer writes are confined to `work/persist_cycle3/`. No repository, manuscript, prior freeze, external account or publication was changed, and no commit was made.
