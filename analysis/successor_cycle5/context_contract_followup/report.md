# Mordstein paired-context contract follow-up

**Decision: NO GO for the proposed outcome branch.** This bounded public-source check recovers a more complete assay-accession map, but not an exact 2020 reporter-CDS index or authenticated independent preparation pairing. The previous reserve decision remains valid. Unknown sequence, matched-unit and independent-preparation counts are null, not zero, in `contract_map.json`.

## What was resolved

The current [ENA study record](https://www.ebi.ac.uk/ena/browser/view/PRJNA596086), secondary accession SRP237851, retains 76 sequencing runs and 32 BioSamples. Their run/sample/library/title mapping exactly matches the prior audit. No processed analysis record was returned by the targeted ENA analysis query. This negative query result applies only to that accession and endpoint, not every possible public deposit.

All eight CM6 experiment XML records were retrieved. Together with four structured BioSamples, they authenticate the following sample categories:

| Context and fraction | BioSample | Library labels | Experiment accessions | Run accessions |
|---|---|---|---|---|
| Intronless, cytoplasmic | SAMN13612104 | CM6_c1, CM6_c2 | SRX7391373, SRX7391374 | SRR10714569, SRR10714568 |
| Intronless, nuclear | SAMN13612106 | CM6_n1, CM6_n2 | SRX7391401, SRX7391402 | SRR10714541, SRR10714540 |
| Intron-containing, cytoplasmic | SAMN13612105 | CM6_ci1, CM6_ci2 | SRX7391375, SRX7391400 | SRR10714567, SRR10714542 |
| Intron-containing, nuclear | SAMN13612107 | CM6_ni1, CM6_ni2 | SRX7391403, SRX7391404 | SRR10714539, SRR10714538 |

The newly checked [SAMN13612106](https://www.ebi.ac.uk/ena/browser/view/SAMN13612106) explicitly identifies HeLa, nuclear fraction, no 5-prime UTR intron and the mRNA-localization assay. [SAMN13612105](https://www.ebi.ac.uk/ena/browser/view/SAMN13612105) identifies HeLa, cytoplasmic fraction and intron presence. This closes the earlier nuclear-intronless sample-attribute gap. All four categories now have direct structured metadata support for HeLa.

The experiment titles still say HEK293, while structured BioSample fields and submitter aliases say HeLa. The conflict is preserved explicitly. No title was silently repaired. The metadata support a same-cell-type assay comparison; they do not establish equal batch, culture history, induction, fractionation or extraction conditions.

## What remains unresolved

1. **Exact CDS units and membership.** The [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7181179/) describes a nominal shared 217-member donor pool and separately states source-category counts summing to 219. It links raw sequencing, not a complete ID-to-CDS reference. The newly inspected deposits do not explain overlaps or exclusions, supply the 2020 membership list, or establish exact variant intersections among assays. No translation or sequence-identity verification is possible from category descriptions alone.
2. **Physical preparation identity.** The eight experiment records describe fractionation followed by GFP-targeted RT-PCR. Their distinct SRX accessions and suffixes 1/2 identify library records. They do not specify independent cell cultures, integration preparations, extraction batches, fractionation pairs, culture dates or randomized matched controls. Shared sample accessions do not resolve whether repeats are technical or biological. All independent-preparation and cross-context pairing fields therefore remain null.
3. **Assay intersections.** No processed ID/coverage matrix was recovered. Nominal library membership cannot substitute for the measured variant intersection after assay filters. No outcome or count table was opened to manufacture this intersection.

The assay map is useful provenance, but is **not** a complete matched reporter/preparation map. Source labels alone do not satisfy the contract.

## Bounded retrieval and failed routes

- The official PMC cloud inventories for PMC7181179 and the cited Mittal 2018 design-origin paper, PMC6112741, were fetched again. Both hashes exactly match the prior inventories. The Mordstein package lists article assets and its known PDF supplements; no reporter FASTA, spreadsheet or analysis code index appears.
- A distinct legitimate fallback succeeded: the [Europe PMC supplementary-files package for PMC6112741](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6112741/supplementaryFiles), 1,397,995 bytes, SHA256 `6bb245ecfede580fa374b8568419ae68ef0c4d51f3cbab77e5fee5c5b6a46fa2`. Only its ZIP directory was inspected. It contains the known supporting PDF and eight figure images, but no `pnas.1810022115.sd01.xlsx`. No image or phenotype table was opened. This does not establish that Dataset S1 is unavailable through every legitimate route.
- The prior publisher Dataset S1 HTTP 403 receipt remains unchanged. That blocked URL was not retried or bypassed. Even recovering the older design file would still require authenticated 2020 membership and assay mapping.
- The installed PMC availability skill used its legacy NCBI endpoint and returned HTTP 404 for both PMCIDs. These are client/endpoint failures, not evidence that either paper is absent. Public cloud and Europe PMC metadata supplied the distinct fallback above. The skill itself was not edited.
- Connected GitHub searches for PRJNA596086 and the Mordstein/codon/splicing repository terms returned no results. The primary data/code statement supplied no code-repository URL. Search absence is not proof that no public author code exists.

Each new request, final URL, status, size, timestamp and saved-byte hash appears in `retrieval_receipts.json` or `experiment_retrieval_receipts.json`. The BioSamples record first-public and last-update dates of 2021-01-01. `prior_hash_check.json` verifies all 31 files in the previous eligibility manifest plus the original paper XML. `build_contract_map.py` independently compares current versus prior run metadata and creates the exact accession map without phenotype access.

## Exposure and stopping boundary

No raw reads, variant-level expression, PSI, read-count outcome table, model or fitted analysis was opened. Ordinary search snippets incidentally exposed already published aggregate conclusions, source-paper correlation summaries and some unrelated cited-paper summaries. They were not used to select a cohort or outcome. The search exposure is recorded in `exposure_ledger.json`; this review is not described as wholly outcome-blind. Targeted local article inspection was restricted to availability sections, supplement titles and previously extracted Methods.

No author contact, protected-access workaround, additional reporter branch, repository modification or publication occurred. The precise inputs that could reopen the gate are an authenticated 2020 variant-ID/full-CDS membership map, assay-specific ID inclusion mapping, and explicit culture/extraction/fractionation preparation identifiers with documented cross-context pairing. Until those exist, the confidence-based paired-context question remains closed. This is an identification and data-contract failure, not a failure of the biological hypothesis.
