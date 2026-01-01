# Pinned human reference proxy, outcomes unopened

Official GENCODE release48 is GRCh38.p14. Its chromosome annotation header identifies GENCODE48/Ensembl114. This is a transparent, fixed public reference proxy. The original CNOT3 reanalysis names GRCh38/Ensembl and a longest-valid-CDS convention but does not identify its release or selected transcript IDs. These files therefore do **not** authenticate the expressed isoform, reproduce the unspecified author annotation, or resolve the original annotation uncertainty.

## Sources and authentication

- Official release page: https://www.gencodegenes.org/human/release_48.html
- Official annotation format: https://www.gencodegenes.org/pages/data_format.html
- Release directory: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/
- Release checksum file: `MD5SUMS`, SHA256 `6c244b5f1672caf7ed5ad49ddc15a0590f1187f2eddfe108e651afa632b85ccc`.
- Cached chromosome-only `gencode.v48.annotation.gtf.gz`:58,971,717bytes; SHA256 `37a298c2b57e41988cfe12667134106523c69c24a1807347679d19e88aa1f4a9`; MD5 `3eb94f7f200b42eb85d56b3afae5b0e7`, matching the release checksum. Existing provenance records the official release-directory URL.
- Downloaded `gencode.v48.pc_transcripts.fa.gz`:52,617,517bytes; SHA256 `d29486ce127cf4ab0f681fdcfa9bcd9fe2cfd18e7331b75257eb12c4cbb385f5`; MD5 `b28c7e02902a96fdd259b3a15ad53014`, matching the same official checksum file. The official page describes this FASTA as coding-transcript nucleotide sequences and lists its included transcript biotypes.
- Whitelisted inventory: `../data_eligibility/reanalysis_gene_ids_only.json`, SHA256 `85fb60552f3c12477d28dbeeafbd171d52828ecb66e30591cefa60f77f482c52`. It contains22,848 distinct unversioned ENSG IDs and cell coordinates only.

Downloads used standard certificate-verified HTTPS, with response metadata saved beside each source. Release MD5 checks establish agreement with the official transport-retrieved manifest; SHA256 pins the retained bytes. They are not a detached cryptographic publisher signature. No assay workbook, HL, IP, quality or conversion-ratio values were opened by this feature workflow. No Ensembl REST call was needed. The applicable Ensembl skill was read before choosing the official release-file route; that skill's REST-script requirement applies if REST calls are made.

## Deterministic selection and coordinate rules

Matching strips only a terminal numeric version from ordinary ENSG identifiers. It does not collapse `_PAR_Y` copies or map aliases. Restrict to chr1–22,X,Y and annotation gene type `protein_coding`. This excludes mitochondrial genes with a different genetic code and does not silently reinterpret noncoding genes or unmapped IDs. Within each eligible gene, consider every transcript supplied in the coding-transcript FASTA, including NMD and LoF biotypes, rather than silently preferring a canonical transcript.

A valid candidate must have an unambiguous A/C/G/T complete CDS of at least6nt, divisible by3, begin with ATG, end with one standard-code stop, and have no earlier standard-code stop. Annotated `cds_start_NF`/`cds_end_NF`, missing complete start/stop features and selenocysteine features are excluded under this declared20-AA/61-sense-codon proxy. No minimum protein-length, expression or assay-value filter is added. Select the greatest CDS length including stop, ties by the literal lexicographic versioned transcript ID. Including versus excluding the one terminal stop gives the same length ranking. All longest-CDS ties are retained in metadata.

GTF exon and CDS/start/stop/UTR genomic intervals are1-based closed. Exons are ordered in transcript5'-to-3' direction, with minus-strand coordinates mapped accordingly. Their projection uses0-based half-open spliced coordinates. Merge CDS and terminal-stop annotations, require one contiguous spliced coding interval, and require complete mapped start/stop codons at its endpoints. Require exact agreement with the independently supplied FASTA `CDS`, `UTR5`, `UTR3` header coordinates and total transcript length. Some GTF UTR features include the stop-codon interval, so subtract the complete CDS including stop before checking the flanking UTRs. UTRs are always from this same selected transcript. Their saved coordinates are spliced coordinates, not genomic coordinates. The original GTF and transcript ID reconstruct their genomic pieces.

The implementation also refuses ambiguous UTR bases; none occurred among the inspected coding FASTA candidates, so that guard did not change the observed longest-valid-CDS selection. Zero-base UTRs have length0 and missing GC, represented by an empty CSV value. `mRNA_start_NF`/`mRNA_end_NF` flags identify separately recorded incomplete UTR annotations even when some bases are present. A numeric GC for such a region describes its annotated bases, not a complete expressed UTR.

## Feature conventions and finite coverage

The feature table has87 numeric columns:61 sense-codon fractions,20 AA fractions, CDS length and GC, and the length and GC for each UTR. Codon fractions divide counts before the terminal stop by the number of sense codons. AA fractions divide the translated protein counts by its length, including the initiating methionine. Both sum to1. CDS length and GC include the terminal stop; UTRs exclude it. GC uses the exact region's nucleotide count. Every selected row has hashes for the transcript, complete CDS, translated protein and both UTR sequences. The table contains identifiers and reference features only, no assay endpoint.

From22,848 target genes,19,811 match chromosome-only annotation. Selection yields19,567 genes. The remaining statuses are3,037 not in this CHR release,142 coding genes without an eligible complete standard CDS,97 annotated noncoding genes and5 non-nuclear genes. These missing IDs have not been remapped or assumed missing at random. All target coding transcripts with GTF CDS features occur in the paired FASTA; zero are missing. The parser stores170,977 target-gene transcript records, of which111,246 occur in the coding FASTA. The59,731 remaining records have no CDS and are not missing coding inputs.

There are81,279 valid coding transcripts. Candidate exclusions are29,519 annotated incomplete CDS,166 lacking complete start/stop annotations,155 non-ATG starts,79 selenocysteine cases and48 internal standard stops. No GTF/FASTA coordinate mismatch or ambiguous sequence was encountered. Selected representatives comprise19,235 `protein_coding`,329 `nonsense_mediated_decay` and3 `protein_coding_LoF` transcripts. The latter two categories are a consequence of the longest-valid-CDS rule and do not imply they are the dominant expressed isoforms. There are4,172 genes with a longest-CDS tie.

UTR5 statuses are18,900 annotated nonempty,662 with no annotated bases and5 annotated incomplete. UTR3 statuses are18,972 annotated nonempty,555 with no annotated bases and40 annotated incomplete. These flags must remain available when root freezes any downstream missingness/coverage rule. No choice about complete-case modeling, imputation, partitions or outcomes was made here.

AA fractions are exact linear aggregates of the61 codon fractions. Thus an AA-versus-codon predictive comparison must account for this nesting and use frozen, training-only fitting/tuning rules. A higher-dimensional in-sample fit alone cannot demonstrate a synonymous effect, biological mechanism or gene-level independence.

## Artifacts, checks and replay

`reference_v1/features.csv.gz` is the reference-only feature table; `gene_inventory.csv` retains every target and its status; `candidate_audit.jsonl.gz` records per-transcript validity; `selected_transcripts.jsonl.gz` retains coordinates, tags, types, sequence hashes and tie records; `summary.json` records complete finite counts and unchanged input hashes. `exact_protein_group_option.json` is an unexecuted partition option, grouping exact translated-protein identities only. `grouping_proposal.md` provides the separately requested near-paralog proposal without aligning any real pair.

Nine synthetic feature tests passed, including independent standard-code agreement for all64 codons, exact fractions, both genomic strands, spliced intervals, incomplete/Sec/ambiguous/invalid CDS, empty UTR semantics, version handling and ties. `verify_selected.py` independently extracts the selected regions from FASTA headers and recomputes translation and all feature arithmetic with Biopython. Its receipt records complete selected-row agreement, sequence hashes and the selection/tie rule. This independently checks arithmetic, not the biological suitability of the chosen proxy. No folding or motif scan is performed.

```text
python -I -B build_features.py --gtf /path/gencode.v48.annotation.gtf.gz --fasta /path/gencode.v48.pc_transcripts.fa.gz --ids /path/reanalysis_gene_ids_only.json --out /path/new-reference-output
```

The three input hashes are fixed in code and output must be fresh. The local independent check expects the named `reference_v1` output next to its sources. The source download script preserves the successful download receipts; an initial local path typo affected only an optional cached-header preview after all three metadata downloads had succeeded. It was corrected before sequence retrieval or feature extraction. No source or scientific output was altered to obtain a match.

This handoff establishes a reproducible reference feature proxy and exposes its coverage and assumptions. It does not authenticate author isoforms, infer causes from half-lives, claim biological family independence or evaluate the CNOT3 predictive hypothesis.
