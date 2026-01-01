# Proposed frozen K562 three-state control audit

Status: **ready for root review and freeze; numeric access remains closed.** No new control outcome, count, RPM, source-score cell or outcome missingness has been inspected. This protocol is prospective only with respect to our extraction and calculation. The source already reports the selected R1207R example qualitatively.

## Question and admissible conclusion

For the published K562 arrayed local states, is combined c.3621A>G plus c.3624G>A readout consistently lower or higher than marker-only c.3624G>A and edit-only c.3621A>G across every admitted source-RPM pairing?

The two contrasts are conditional-background comparisons of source-defined local-word cDNA/gDNA ratios: combined/marker-only and combined/edit-only. They are not WT-isolated intended-edit effects, interaction contrasts, complete haplotype effects, population-average biological effects or independent replication. WT absence does not make these narrower contrasts undefined. It does prevent the missing factorial inference.

## Fixed source and metadata contract

- Source: official Supplementary Data 6 of DOI 10.1038/s41467-026-78344-9, extracted from the already authenticated official Supplementary Data 1-9 ZIP.
- Workbook: work/research_cycle8/cprime_gate/supp_Supplementary Data 6.xlsx.
- SHA256: **a3d18aff825b756711de622533cb730c1202e571f86a3b9b2e00edea6b073dd3**.
- Official ZIP: https://pmc-oa-opendata.s3.amazonaws.com/PMC13675207.1/41467_2026_78344_MOESM3_ESM.zip.
- ZIP SHA256: **a88a9850874d0ef21478406c6e7c9e437dc7fe7f29ece32c84033efe8940a614**.
- Worksheet: **validation K562**. It contains all three selected labels. Other worksheets and all other genotypes are outside the outcome read whitelist.

| State | Exact source label | Row / label cell | gDNA RPM cells | cDNA RPM cells |
| --- | --- | --- | --- | --- |
| Marker-only M | R1207R_ATT(G/A) | 15 / C15 | F15,G15,H15 | I15,J15,K15 |
| Edit-only E | R1207R_(A/G)TTG | 16 / C16 | F16,G16,H16 | I16,J16,K16 |
| Combined A C | R1207R_(A/G)TT(G/A) | 17 / C17 | F17,G17,H17 | I17,J17,K17 |

F/G/H are labeled K562_gDNA_1_RPM, K562_gDNA_2_RPM, K562_gDNA_3_RPM. I/J/K are labeled K562_cDNA_1_RPM, K562_cDNA_2_RPM, K562_cDNA_3_RPM. Exactly these **18 numeric cells** may be read after root freezes and authorizes the run. FC columns L-N, all other numeric cells and raw reads remain closed.

| State | Intended local substitutions | Identical deposited DNA and cDNA search word |
| --- | --- | --- |
| M | c.3624G>A | GGTTTGAGAGATTAGTGGAACTGG |
| E | c.3621A>G | GGTTTGAGAGGTTGGTGGAACTGG |
| C | c.3621A>G; c.3624G>A | GGTTTGAGAGGTTAGTGGAACTGG |

The common local WT reference word is GGTTTGAGAGATTGGTGGAACTGG. Its changed positions are 11 for E, 14 for M and both for C. The annotation is ENST00000358273; these are reference-synonymous substitutions. Local words and source annotations do not establish complete realized genomic haplotypes or the identity of every processed product.

The metadata contract lists all 18 matching ENA runs, their BioSample accessions, modality and source suffix 1/2/3. The selected condition-name families are:

- E: K562_wo_marker followed by R1207R_(A/G)TTG and modality/replicate.
- M: K562_wo_edit followed by modality and R1207R_ATT(G/A)/replicate.
- C: K562_marker followed by R1207R_(A/G)TT(G/A) and modality/replicate.

There are six labeled runs per state. **Twelve earlier unprefixed technical-family runs are excluded as a distinct preparation family**, not pooled or averaged into these observations. No KBM7 outcome is admitted. GtoT is excluded because its deposited counted motif equals the A-marker combination; GtoC is outside this fixed three-state comparison. The contract enumerates exclusions before numeric access.

## Preparation meaning and eligibility restrictions

The source describes arrayed K562_PE_MLH1 cells and triplicate frozen aliquots, with 31/21-day harvest wording. It does not resolve harvest day, clone and independent transduction/culture identity for each arm. Named ENA condition/modality/suffixes are compatible with the three source labels, but no explicit source-table-cell to BioSample/run join was recovered.

Consequently, source-index matching is an operational convention, not authenticated physical pairing across arms. Three columns are not assumed to be three independent cultures. Exact condition-name compatibility does not supply a batch-randomized or paired biological comparison. The permitted question is finite-source descriptiveness conditional on the source's common arrayed-assay interpretation. Any new sample mapping or contradictory source record must be reported before claiming biological comparability.

## Extraction and arithmetic

1. Refuse workbook/hash, worksheet, label-cell or header mismatch before numeric extraction. Refuse optimized Python execution and existing output paths. Preserve every failure.
2. Read only the 18 allowlisted numeric cells from Excel XML. Convert stored decimal numeric strings to exact rationals; do not use display-rounded values or silently change the source scale. If a cell is a formula, only a present stored numeric cache is admissible; do not execute Excel or infer a missing cache. Record the lexical source and cell coordinate.
3. Do not renormalize within the selected three states, add a pseudocount, recalculate source RPM totals, trim values or replace missing cells. Use deposited RPM as deposited. Its normalization and unobserved recovery remain part of the source measurement process.
4. For each state s and indices i,j in {1,2,3}, form q(s,i,j)=RNA_RPM(s,i)/DNA_RPM(s,j). Retain all **nine** constructed possibilities per state and their extrema. Also list the three source-index ratios i=j, clearly as an operational diagnostic, without narrowing the primary envelope.
5. For each comparator S in {M,E}, form all **81** quotients q(C,i,j)/q(S,k,l). Retain every case and the exact lower/upper envelope. The simultaneous two-contrast completion space has **729** triples of selected state ratios, sharing the same C value. These are sensitivity completions, not 81 or 729 independent observations.
6. No means, medians, fitted model, p-value, confidence interval, hypothesis test or selected replicate subset enters the primary analysis. The exact identity of extrema is reported, with all cases available for review.

## Missingness, zero and schema failures

- Missing, nonnumeric, uncached formula, nonfinite or negative RPM: preserve the offending cell and return a source/schema or unknown-measurement status; no silent exclusion or replacement. A state with an unknown entry has an unbounded [0,+infinity] primary envelope and affected contrasts are unknown.
- A zero DNA RPM makes the affected state ratio undefined. Retain the state and make its primary envelope unbounded; do not repair it with a pseudocount.
- Zero RNA RPM with positive DNA gives a recorded-data ratio of zero. This is not evidence of a biological zero.
- In a contrast quotient, positive/zero is recorded as +infinity with an assay-zero flag. Zero/zero is undefined and yields an unknown [0,+infinity] contrast envelope. Do not convert undefined cases into either direction or a finite fold estimate.
- If every primary entry is defined, exact rational and extended-nonnegative comparisons determine the envelope. Infinite endpoints and every zero/unknown flag remain visible. No reported sign is a calibrated biological sign.

## Frozen comparisons and complete reporting

For each C/M and C/E envelope [L,U], report all of:

1. Exact L and U, with finite, infinite or unknown status.
2. Strict directional separation at one: uniformly lower if U<1; uniformly higher if L>1; otherwise neither strict direction. Equality/touch at one is not strict separation.
3. Secondary at-least-twofold display only: uniformly at least twofold lower if U<=1/2; uniformly at least twofold higher if L>=2. Exact boundary equality is separately flagged. This does not replace the primary at-one statement or introduce a biological practical tolerance.
4. The source-index diagnostic and whether it agrees with the full Cartesian result.
5. All three state envelopes, all 27 state-ratio possibilities, all 162 contrast possibilities and any unknown, zero, disagreement or null. No state, comparator, cell family or sensitivity is chosen after outcomes.

The primary strong descriptive conclusion requires the corresponding strict inequality for **every** admitted completion, not existence of a favorable assignment. If only one conditional contrast separates, report that asymmetry. If neither does, retain the informative failure. A successful C/M separation would reject equality of those two deposited local-state readouts under this finite pairing envelope. It would not prove a pure intended-edit effect, because the intended edit is compared in the marker background and source recovery/background can differ. C/E has the corresponding marker-background qualification.

## Prior, exposure and validation

The original paper already reports selected R1207R combined-state behavior and effects of marker removal, with Supplementary Fig. 1e/f and Data 6. This audit cannot discover marker non-neutrality, claim an untouched biological validation set or independently confirm a molecular mechanism. The potential increment is an exact, complete-range check of two conditional-background comparisons after excluding the unresolved sequence alias.

Qualitative source conclusions, source selection of this example and metadata were exposed before this protocol. **The 18 proposed numeric cells and their missingness remain unparsed by this Cycle9 audit.** Numeric access requires root's protocol/code freeze and explicit authorization. An independently written arithmetic checker should verify every admitted quotient and boundary after that run; this is numerical verification, not biological replication.

Stop after one unchanged calculation and its independent check. Do not add a new threshold, move to another example, pool the earlier technical family or access another source to rescue a weak contrast. Root decides whether the likely narrow contribution justifies opening these already source-selected controls.
