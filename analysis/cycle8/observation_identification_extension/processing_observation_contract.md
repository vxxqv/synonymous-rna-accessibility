# Dao: what the documented observation architecture identifies

**Empirical decision: no authenticated full-transcript-to-coding-product matrix is available in the inspected records.** The architecture does identify operational read-category contrasts, conditional on its filters and assignments. It does not identify the true coding-product mixture or make the large regulatory library a synonymous coding landscape. The small matrix below is a declared abstract refinement of the observed categories, not a fitted matrix or an assertion that its latent species occur in the source.

## Authenticated rows and unknown columns

The originating journal article is Dao et al. 2025, DOI 10.1038/s41467-025-62000-9, PMC12297286. The author repository is MustoeLab/Publications at 7a8779a484af797953a154966e8b06b4afced3ed. Existing reference auditing establishes one common 720 nt intended ORF, while the variable reporter elements occur downstream. Reference protein identity does not resolve altered processed products.

The source Methods describe RT using a downstream primer and amplification between a GFP-side primer and a downstream element-side primer. A required forward-primer match and barcode assignment precede alignment/category construction. The readout therefore conditions on recoverable amplicons and retained reporter identities. Missing primer/barcode classes, artifact categories and high-mismatch classes are not a census of every transcript. The paper also describes protocol-dependent suppression of palindrome-containing templates. Equal recovery of all RNA products cannot be inferred from shared primers or a DNA library control.

| Layer | Authenticated observation definition | Unknown information |
| --- | --- | --- |
| Raw per-preparation category counts | Counts by assigned reporter/barcode and `best_category`; accepted full-length class and individual low-mismatch reported junction classes are retained. | Actual clone/genotype origin for every assigned read, undetected RNA products, misassignment/capture probabilities and complete transcript path. |
| Junction classifier | `categorize_read` sums recognized GT-AG deletion lengths for mismatch adjustment but keeps donor/acceptor coordinates from the last recognized deletion in its loop. The combined category retains that reported junction label. | Earlier junctions are not encoded in that label. Full CIGAR/path can carry more information than the final category. A reported junction label is not by definition a full isoform. |
| Normalized raw fractions | Within each reporter in one preparation, eligible category count divided by the eligible count sum. | Absolute RNA abundance; category-specific RT/PCR recovery; excluded/unobserved molecules. |
| Deposited median fractions | Barcode/category pivot after preparation concatenation, followed by a coordinatewise median over internal reporter barcodes. | A normalized isoform-probability vector, independent biological replicates, or a linear observation equation for molecular amounts. |
| Coding product | Known intended reference translation and source-reported mechanisms capable of changing processing. | A complete clone/barcode/full-junction-path/CDS/stop assignment for each observed molecule. No translated-product weights can be filled in from the category name alone. |

Raw counts, source preparation identities and some read-level files are available by source metadata, but no such numerical observation was read or analyzed here. The current task does not infer actual category presence, missingness or a per-source matrix rank from a Methods schema.

## The identified contrast is conditional and operational

If retained raw counts for one assigned reporter/preparation are F and S_j for its reported junction labels, then

    q_j = S_j / (F + sum_j S_j)

and the corresponding full-length fraction, sums, odds or differences between those conditional read fractions are determined by those counts when the denominator is positive. This remains a useful read-defined RNA-processing contrast despite unknown full clone/product mapping. It is not automatically the fraction of all original RNA molecules spliced at that site, a transcript-specific export rate or a protein fraction. Preparation contrasts inherit the source selection and calibration; counts and barcode labels do not supply independent molecular or biological units by themselves.

Even in a fixed linear expectation idealization, physical transcript amounts x_t would produce expected reads

    y_j = sum_t q_t * 1{reported_category(t)=j} * x_t,

where recovery factors q_t and category assignments must be known to specify the matrix. Mapping errors would require a more general assignment kernel. Those coefficients are not authenticated here. PCR need not satisfy that idealization. We therefore do not give this formula an empirical numerical rank.

## A small conditional matrix and explicit ambiguity

At one assigned reporter, consider one full-length category F and one SPECIFIC reported junction bucket J, not a merger of all junctions. Suppose its unresolved latent refinement contains a full-length species, two distinct paths with the same reported J and an unretained species. This possibility is not excluded by the missing path/clone map. As captured-read masses, their aggregation is

              F  J_alpha  J_beta  U
    A = [     1     0       0    0 ]
        [     0     1       1    0 ].

This abstract matrix has rank 2 and kernel basis (0,1,-1,0),(0,0,0,1). Every other category row can be held unchanged; the ambiguity is within one reported junction bucket and is not manufactured by discarding additional observed junction labels. It says nothing about how many hidden species actually occur in the source. Under the additional favorable idealization of known equal unit recovery for retained transcript species, the same A acts on hypothetical transcript abundances; U is an excluded species with zero recovery. This idealization is not empirically authenticated and is used only to show ambiguity before additional unknown recovery is introduced.

Two strictly positive rational mixtures are

    u=(2,3,1,1), v=(2,1,3,1).

They both have total 7 and observation vector(F,J)=(2,4). If the first path and full-length species preserve one coding product while the second path/unretained species do not, coding contrast l=(1,1,0,0) takes values 5 and3. Even an extra total-amount row does not resolve this null direction because(0,1,-1,0) has zero sum. The coding labels are a hypothetical unresolved assignment, NOT authenticated source isoforms or a claim that such mixtures were observed. If the actual product assignment is constant within every category, this particular ambiguity disappears; recovering that fact is precisely the missing mapping obligation.

The relevant general linear statement is elementary: for a known A, a linear contrast l*x is globally determined by y=A*x on an open positive domain if and only if l lies in the row space of A. Sufficiency follows by writing l=b*A. For necessity, l outside that row space has nonzero pairing with some kernel vector; a sufficiently small positive and negative perturbation along it preserves a strictly positive x and y but changes the contrast. At a boundary-specific observation a contrast can happen to become fixed; that does not establish the global statement or authenticate A. The actual source's A and coding contrast are unknown, so this lemma supplies a missing-information test, not an empirical identification result.

## Source-parser evidence for category compression

An isolated AST extraction of the original `categorize_read` function was tested on symbolic alignment predicates. CIGARs `100M30D50M` and `10M30D60M30D50M`, with their matching deletion edit distances and stipulated valid deletion-boundary predicates, both return the same final donor/acceptor label at 100/130. The second CIGAR has an additional earlier deletion. No nucleotide sequence was generated, and no claim is made that these two paths or boundaries occur in the actual reference/data. The test demonstrates loss of earlier-path information in this software output, not an observed biological deletion or a new molecular mechanism.

This is narrower than alleging a source bug. The reported category can serve the source's chosen endpoint. It is simply insufficient by itself for a new full-path/coding-product interpretation. Full read alignments may retain information omitted from the count table and would require their own approved mapping audit, not a presumption that the data are irretrievably lost.

## Why the processed median table cannot be treated as A*x

Three source-code details matter before any physical interpretation:

1. `process_raw_count` retains full-length/spliced categories, rejects per-reporter totals<=10, and normalizes within the retained reporter count. Excluded counts are not included in that denominator.
2. `filter_barcode_count` uses the number of full-length rows per RE after concatenating preparations, rather than explicitly counting distinct barcode identities. Thus three barcodes observed in two preparations can supply six rows and pass a threshold 5. This is an implementation/counting fact, not a claim about the actual retained source population.
3. `pivot_table` defaults to averaging repeated barcode/category entries across preparations, then the source computes separate category medians over reporter barcodes. Missing categories are filled after pivoting, not explicitly zero-completed per preparation first.

A synthetic eligible example uses six barcodes with 100 reads each and full-length counts present. Two barcodes each have category fractions(0.8,0.1,0.1), two(0.1,0.8,0.1), and two(0.1,0.1,0.8). The three coordinatewise medians are each 0.1 and sum to 0.3. Therefore even perfectly normalized per-barcode input vectors need not produce a probability vector after the median step. In a separate two-preparation example, fractions(F,S)=(0.9,0.1) and(1,0), with the second S row absent, give pivot values 0.95 and0.1 rather than the zero-completed preparation mean 0.95 and0.05. These are algebraic tests of the documented operator, not empirical artifact frequencies or an accusation that its intended statistic is invalid.

Neither operation can silently be replaced by a pooled count ratio, a biological-replicate likelihood or a normalized transcript mixture. An intended per-barcode/per-preparation read-category statistic remains definable, but needs its exact operator and dependence stated. No alternative endpoint is selected here.

## Biological decision and missing record

Retain the source-defined read-category contrast as an operational observation. Close the stronger empirical coding-product identification claim at the metadata gate. The most informative missing existing record is a linked reporter/barcode-to-realized-full-plasmid-to-complete-RNA-junction-path/CDS/stop assignment, with preparation identity and capture exclusions retained. If the question is a physical RNA fraction, relative recovery/calibration information is additionally required. The known intended FASTA, a category label, a barcode consensus or a selected Western blot example does not supply that record.

This is a source/observation limitation, not proof that the source cannot ever support a processing analysis. The paper already establishes cryptic splicing and product changes; no rediscovery or new synonymous-neighborhood claim follows. No phenotype count, outcome table, raw read, new sequence, experiment, intervention or clinical inference was introduced.
