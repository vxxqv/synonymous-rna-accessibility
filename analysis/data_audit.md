# Bounded source and coordinate audit

## Mustoe 2018 processed mirror

The study used E. coli K-12 MG1655 and aligned to GenBank U00096.2, as stated in the published methods. Current U00096.3/NC_000913.3 coordinates must not be substituted silently. All 194 full transcript sequences in the pinned ShaKer mirror map uniquely and exactly to U00096.2: 100 plus strand and 94 minus strand. Each condition contains 442,421 nucleotides. The same full sequence is present in all three conditions.

The audited mirror is [BackofenLab/ShaKer at commit 6f29547d70a2b72ceccf28218c790c45406946e8](https://github.com/BackofenLab/ShaKer/tree/6f29547d70a2b72ceccf28218c790c45406946e8/data/weeks194_orig). The full transcript DBN/reactivity files, conversion script, and checksums were retrieved. Every DBN is balanced and every decoded pair is one of AU, UA, GC, CG, GU, or UG. Decoded pair counts are 120,830 in-cell, 130,540 cell-free, and 132,193 kasugamycin. Missing reactivities number 51,235, 67,034, and 192,510 respectively.

The mirror converter converts a CT partner into an opening or closing parenthesis based solely on endpoint order. Crossing CT partners cannot generally be recovered by stack-decoding that representation. The original 2018 archive remains unavailable: the historical laboratory archive returned HTTP 403; the document-linked webshare archive had an expired TLS certificate and returned HTTP 404 after the certificate issue was diagnosed. The published methods describe SuperFold/RNAstructure 5.8 and no pseudoknot insertion, but SuperFold supports manual insertion, so that description alone does not establish byte-level or partner-level identity. The mirror is therefore an explicitly derived noncrossing structural model until original CT correspondence is verified.

The `*_genes.dbn` files are cropped sequences refolded by the mirror's `remove_genes.py`; they are not the original published transcript pair tables. That script also contains a start-coordinate slicing ambiguity. The analysis uses only the full transcript records and independently mapped annotation.

## Coding map and eligibility

The frozen U00096.2 GenBank annotation contains 4,321 CDS features, including complex joined or interrupted CDSs. The simple contiguous-feature parser maps 409 fully contained same-strand CDS occurrences, representing 407 loci. Complex CDSs are outside the one-frame model. Every mapped simple CDS length is divisible by three and has a terminal stop. Standard internal translation agrees with the annotated protein except fdhF/b3894, which uses selenocysteine recoding and is excluded from the standard-code analysis.

Initiators include 367 AUG, 31 GUG, 10 UUG, and one AUU occurrences. Freeze the native initiating codon rather than treating it as an ordinary internal synonymous codon. Freeze terminal stops. Seventy-seven mapped CDS occurrences overlap another annotated CDS; variable blocks must avoid affected genomic positions, including opposite-strand obligations.

Mapping preserves duplicated loci: rpsO/b3165 appears in two overlapping transcript contexts, as does acpP/b1094. Window-level results and duplicated contexts must not be treated as independent biological replicates. Transcript/operon dependence and genomic overlap remain inferential considerations.

Finite negative SHAPE reactivities are valid background-subtracted measurements. Coverage uses missing versus finite values, not a sign cutoff. Missing reactivity does not establish missing reads or zero structure. The mirror lacks original per-position read-depth and uncertainty columns, so finite-value coverage cannot recreate the original study's greater-than-1,000-read eligibility criterion. The independent original-CT branch instead preserves full ShapeMapper profiles.

## Independent original-CT branch

The CC0 PAIR-MaP 2019 source provides original CT files, direct principal/minor interaction windows, and full probing profiles. Its audited rpsM/rpsB branch is documented in `analysis/pairmap/README.md`, with 27 selected original files, exact sequence mapping, a one-to-many coding-position ledger, and 40 complete condition/window graph results. Across 150 nonempty filtered graphs there is no additional fragmentation. Ten filtered sets are empty and their accessibility metrics remain undefined. No qualifying one-to-three-nucleotide coding register-shift example was found. These negative results are retained explicitly.

## Reporter exclusion decision

A pinned processed Goodman 2013 workbook was recovered from [SalisLabCode commit 46da9132740bce30bccd107f02808adf2215e2f0](https://github.com/hsalis/SalisLabCode/blob/46da9132740bce30bccd107f02808adf2215e2f0/ModelTestSystem/examples/RBS/datasets/Goodman_Science_2013.xls). SHA-256: `1295f0761b70afc6f6a23c1c3c674f011a1cf156728d850ec04a3c21f65571b9`. It contains exactly 14,234 constructs from 137 genes and 1,096 promoter/RBS/peptide backgrounds; 1,082 groups have 13 variants and 14 groups have 12. RBS sequence and N-terminal peptide are constant within each declared background. Published protein range flags retain 7,327 constructs; requiring sufficient DNA and RNA counts retains 7,316.

The mirror's context-reconstruction code incorrectly treats the continuous `Offset.RBS` column as a transcription-start coordinate. The original publisher table was inaccessible, and a complete transcript boundary was not independently validated. Explicit truncated contexts can be constructed, but a confident external prediction validation would require resolving their biological boundary. Reporter model fitting and held-out prediction were therefore not performed. No incremental prediction result is claimed. The workbook is retained only in the working audit area pending source-license review, with its URL/checksum audit retained in the repository.

For Kudla 2009, a pinned secondary R workspace yielded 154 fluorescence/feature rows but no corresponding sequences. Its source table attributes those data to personal communication. The sequence-to-phenotype gate remains unresolved; it is not a substitute for genotype-level validation.
