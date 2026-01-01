# Cycle9 cPRIME control recovery and narrower comparison

**No new matched WT RNA control or corrected T-marker identity was recovered in the bounded official-source audit. A qualified K562 three-state descriptive comparison is now metadata-resolved and ready for protocol review, with numeric access still closed.** It can compare combined versus each single local state; it cannot supply a missing WT-isolated effect or interaction.

## What was actually checked

1. Reused the earlier source audit, exact sequence metadata, official supplement descriptions, primary Methods and 330-row ENA run-metadata inventory. No earlier failed publisher, PMC-challenge or supplementary timeout route was retried.
2. Queried the current author repository, exact default-branch HEAD/tree and the linked Zenodo latest-version endpoint. HEAD remains 30125415304664c33370cd7801f820a345834e40. All **74 paths, including 63 blob files**, have identical identities to the prior inventory. The current actual tree SHA is 5135ee38a7f24a5e5ae3dc844c10f1a8af0fb4d6; the earlier response's commitish SHA is not evidence of changed contents. Zenodo still returns record 21991157 and the same v1.0.0 ZIP, MD5 c2a125f3d630703032100a8bd548b948.
3. Used the connected GitHub search for the odds-ratio source. It returned no matches; that search absence was not treated as file absence. Ordinary pinned primary GitHub retrieval succeeded for Odds_Ratio.R and all four cDNA/gDNA count-table containers.
4. Used Scholar Sidekick's current correction/retraction check. It reports no correction, concern or retraction on record for DOI 10.1038/s41467-026-78344-9. A targeted primary-domain search identifies the same publication and supplements. This does not guarantee no undocumented correction anywhere.
5. Inspected only first-row schema, second-row identity headers and column-A Sample IDs in the newly acquired control workbooks, plus source code. Their numeric observations and missingness remain closed.

The source/access ledger gives exact URLs, HTTP status, SHA256 and, for the five newly acquired GitHub files, verified Git blob SHA1. All new ordinary requests succeeded. Historical access limitations are recorded separately rather than relabeled as absent data.

## Findings that matter

| ID | Finding | Evidence and interpretation |
| --- | --- | --- |
| C9-F1 | No new version in the checked author-linked routes | Current default branch and latest Zenodo metadata are unchanged. This is a bounded check, not a claim that no other public material could exist. |
| C9-F2 | The apparent WT control is gDNA-only in these deposited odds-ratio containers | Both gDNA workbooks have six control sheets, including WT in edited/unedited cells. Each cDNA workbook has one edited-assay sheet with 98 variant identity rows and no separately labeled WT cDNA state. Odds_Ratio.R loads only the gDNA sheets for its edited/unedited signal-to-background comparison. It is not a WT cDNA/gDNA response. |
| C9-F3 | The source labels “without marker” in this background analysis do not establish separately installed edit-only RNA samples | Those sheets count alternate words in genomic material. The separate arrayed Data 6 edit-only conditions are the relevant RNA controls. Do not substitute source gDNA editing-signal odds for a transcript response. |
| C9-F4 | The unresolved GtoT count-reference mismatch persists | Its deposited DNA and RNA word remains identical to the original G>A combined state. The repository is unchanged. This is a count-reference/label discrepancy, not proof about the physical edited cells. Exclude it; do not repair the sequence from its phenotype. |
| C9-F5 | Three unambiguous local states are represented in K562 Data 6 | Worksheet validation K562 contains marker-only at row 15, edit-only at row 16 and combined G>A at row 17. Source labels and local words are matched in three_state_contract.json. These are exact source identities, not complete realized haplotypes. |
| C9-F6 | Condition-name provenance is usable but biological pairing remains incomplete | The later K562 named families have six RNA/gDNA runs per state, source suffixes 1-3. Twelve earlier unprefixed technical-family runs are distinct and excluded. The source does not provide the arm-level culture/clone/harvest/RNA-DNA-table join needed for a physically paired biological comparison. |
| C9-F7 | A narrow additional calculation is identifiable on the deposited scale | Conditional source C/M and C/E ratio envelopes do not require WT. Full Cartesian RNA/DNA pairing protects against assuming the undocumented pairing. It does not protect against unknown biological/background/recovery differences or establish independent replication. |

The complete 330-row cached ENA metadata has no WT/unedited cDNA-labeled run in the control-name search. That supports the scoped missing-record verdict, not an assertion that unedited RNA could not be recovered from raw reads. Raw reads and outcome fields were not opened.

## Exact eligible object and restrictions

The fixed candidate is **only the three K562 rows15-17** of official Supplementary Data 6, workbook SHA256 a3d18aff825b756711de622533cb730c1202e571f86a3b9b2e00edea6b073dd3. DNA/RNA local words differ from the common reference at positions11 and/or14, encoding intended c.3621A>G and marker c.3624G>A. Their identical deposited DNA/cDNA search words resolve the relevant local identity, while distant edits, complete haplotypes and all processed products remain unauthenticated.

The paper describes K562_PE_MLH1 for the arrayed assay and freezing aliquots in triplicate, with 31/21-day wording. This does not establish three independent cultures or transductions. Matching RNA/gDNA suffixes and ENA condition strings are source conventions, not evidence of common physical aliquots across different states. The comparison is therefore descriptive of deposited local-word readouts conditional on the source's common assay interpretation.

K562 is preferable to KBM7 for this bounded comparison because all three labels appear together in the K562 validation worksheet and its later named families. The KBM7 single/combined records span technical-validation and validation tables without a resolved batch join. No cell lines or earlier/later K562 families are pooled.

## Closest prior and possible increment

The [original cPRIME paper](https://doi.org/10.1038/s41467-026-78344-9) already reports that selected marker combinations alter the RNA readout and compares marker removal in the R1207R example. The [author's source](https://github.com/Boettcher-Lab/cPRIME/tree/30125415304664c33370cd7801f820a345834e40) supplies its count-reference words and signal-to-background analysis; the [linked archive](https://doi.org/10.5281/zenodo.21991157) is unchanged. Replotting the example or finding a marker difference is not a new biological discovery.

A distinct but limited question is whether **both** combined-versus-single comparisons strictly separate across every admitted pairing of the deposited RPM observations. This tests robustness of two conditional-background contrasts, not interaction or intended-edit invariance. The frozen proposal retains all three states, both directions at 1, the unchanged twofold secondary display, all unknowns and the complete null/failure outcomes. It uses no new model or outcome-selected subset.

## Gate and smallest next action

- **Qualified GO:** root may freeze and independently review three_state_protocol.md and its metadata contract for the finite-source K562 comparison. It admits exactly18 numeric RPM cells after explicit authorization, all 9 within-state RNA/DNA pairings and all 81 possibilities per conditional contrast. No numeric cell was read to choose this rule.
- **HOLD:** no causal intended-edit, interaction, biological confidence, independent-validation or full-haplotype claim. No matched WT RNA response was recovered.
- **Precise missing record:** a source specimen manifest linking each of these three arms to clone/background, transduction/culture batch, harvest day, RNA/gDNA extractions, Data 6 cells and run accessions. A compatible WT RNA/gDNA record is additionally required for the absent factorial comparison; it is not required for the narrower descriptive comparison.

If root proceeds, freeze the selector, exact source hash, endpoint, zero/missingness rules, complete pairing space and numerical checker before opening the 18 cells. Report the resulting complete calculation whether it supports one, both or neither conditional contrast. Do not search another example or change tolerances to improve it.

## Scope/object filters and exposure

Scope filter: no broad source search, raw-read reprocessing, new constructs, clinical interpretation or new biological model. Object filter: retain source-defined local words; exclude the T-marker alias, alternative C-marker state, other cell line and earlier technical-family duplicates. No matched WT is inferred from the reference sequence, gDNA background, total read count or a nominal cPRIME value of 1.

New reads were metadata, source code and primary prose only. The primary prose's qualitative selected-example conclusions and published prediction/classification aggregates were exposed; they are not untouched validation. The source numeric RPM cells, source FC cells, their missingness, notebook outputs and raw reads remain closed. The initial combined console read was truncated; the user-request tail and methodology documents were subsequently read separately. One metadata-print helper raised KeyError on a directory entry lacking size; it did not inspect outcomes, and the successful comparison uses path/type/blob identity.

Research progress: canon, charter, bounded source collection, scope filter, object filter and source-backed conclusions are complete. This audit ends before outcome access.
