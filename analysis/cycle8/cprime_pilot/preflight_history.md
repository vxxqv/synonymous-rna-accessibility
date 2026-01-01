# Pre-outcome implementation corrections

No deposited phenotype or count value was read during these corrections.

1. The metadata-only selector initially attempted substring membership on a null consequence annotation. It stopped before creating eligibility. Null annotations now remain outside the annotated synonymous set. The resulting30 groups and66 states independently match B0. This changes no phenotype or scientific threshold.
2. Metadata show the30 eligible multi-marker groups occupy five exons (17,25,27,39,42), not all six source exons. The draft's six-exon wording was corrected before the freeze. Exon34 lacks an eligible multi-marker synonymous group under the fixed annotation rule; it was not removed on outcomes.
3. The first synthetic selective-reader run completed its numerical checks but Windows rejected deletion of its temporary workbook because a ZIP member handle remained open. The reader now parses an in-memory copy of the sheet XML, closing the archive handle. The failed temporary fixture is retained locally. The successful rerun records1,100 exact range comparisons and500 independent shared-order permutation comparisons. No source workbook values were accessed by either run.

These are development corrections, not amendments to an already frozen scientific test. Independent review and an immutable freeze must precede production outcome access.
