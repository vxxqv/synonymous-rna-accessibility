# Exact S4 recovery and coding-gate effect

**Recovered; coding verdict unchanged.** The earlier Springer CDN DNS failure did not indicate a missing or restricted supplement. The official PMC Open Data inventory lists the exact S4 workbook, and one ordinary public retrieval returned the checksum-matching file. It contains promoter/upstream-intron sequences, not an annotated full reporter-CDS/isoform map. This closes the specific recovery question but does not authenticate same-protein synonymy.

## Provenance and bounded route

The already inspected primary XML for [Shanas2026](https://doi.org/10.1038/s41467-026-76370-1) identifies `41467_2026_76370_MOESM6_ESM.xlsx` as Supplementary Data4; the main text references it for promoter constructs. A request using the installed PMC skill's older OA-service URL returned404. This was recorded, not retried or treated as permission denial.

The normal public [official PMC inventory](https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix=PMC13534650&max-keys=100) returned200 and listed the exact object. Only that identified workbook was downloaded, using [its listed public link](https://pmc-oa-opendata.s3.amazonaws.com/PMC13534650.1/41467_2026_76370_MOESM6_ESM.xlsx). No other supplement, figure, outcome table or archive was fetched.

- Exact file size: 11,971 bytes, matching the inventory.
- MD5: `29902e2f98b9db0ef88781246ee32aa7`, matching its single-object ETag.
- SHA256: `9acf21d4b939451248dbebb911fc3769228a079ea64da6404db49779eeac316e`.
- Workbook: one sheet, `Sheet1`, with only `promoter construct` and `sequence` headers.

The existing selective spreadsheet reader first inspected sheet names/headers. After the metadata-only title/schema was confirmed, only those two columns were decoded. No new author code ran. The selected metadata contain 11 entries: pEF1a, pACTB, full/short/minimal CMV, two full-CMV upstream-intron versions, two MYC promoters and two TP73 promoters. Sequence lengths range from 46 to 1,703 nt.

The sequence schema's `non_DNA_characters` field means characters outside **ACGTN** after uppercasing; it does not assert unambiguous bases. pACTB contains one ambiguous N, retained explicitly in the corrected schema summary. No base was repaired or inferred. The preceding summary snapshot is retained as history; source sequence bytes did not change.

## Exact effect on the gate

These are promoter and upstream-intron metadata. They do not annotate a complete shared GFP/cassette coding sequence, coding phase, start/termination boundary or included/skipped isoform CDS. Absence of the exact GFP anchor used for a schema check is not, by itself, a general proof about every possible GFP sequence; the source table's construct scope and lack of coding annotation are the relevant limitation.

Consequently:

- The prior **file-access** gap is resolved through an ordinary official public source.
- The **full coding/isoform identity** gap remains. Identical cassette IDs still do not prove identical translated products.
- The earlier all-2-modulo3 cassette-length finding and conditional coding argument are unchanged. Synonymous graph counts remain “not identifiable / not computed,” not zero.
- Preparation lineage, actual measured context overlap and a distinct endpoint remain governed by the completed original contract report. No outcome values or availability were accessed to revise them.

The original `report.md`, `evidence.json`, selector and their frozen receipts were not modified. This follow-up is additive only. No further retrieval or phenotype analysis is recommended as part of this bounded task; root can close the same-protein route while retaining conditional, explicitly re-scoped matched-cassette portability as a separate research decision.
