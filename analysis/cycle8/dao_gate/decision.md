# Dao processing gate: authenticated object and current limit

The journal version is Dao et al., Nature Communications 16,6844 (2025),
DOI10.1038/s41467-025-62000-9, PMID40715118, PMC12297286. The older
PMC11326173 record is a preprint and is not the authoritative version.
The publisher page failed through the web reader; the DOI-matched PMC
paper and author repository were inspected. The paper links PRJNA1116243
and Zenodo15700943. Author code/data are stated to be MIT licensed.

The pinned author repository is MustoeLab/Publications,
commit7a8779a484af797953a154966e8b06b4afced3ed. Its complete tree, relevant
READMEs, junction-parser source and processing code are retained. Notebook
outputs and individual phenotype rows were not read. Raw file acquisition
preserves exact Git blob hashes; the connector could not retrieve the
binary FASTA, so a direct public author-file download was verified against
the tree blob instead.

## Concrete source result

The intended-transcript FASTA has 6,500 records of length1,045nt.
All contain the same first720nt open reading frame, encoding239 amino
acids plus stop. The first varying nucleotide is zero-based index741,
after that reference coding sequence. Grouping only the final numeric
barcode suffix gives642 reporter labels:640 have10 barcodes, whereas
BBBB and7pcx4 have50 each. These are sequence/index counts, not counts of
independent biological units or usable phenotypes.

This library therefore supplies a regulatory reporter object, not a
library of changing synonymous coding nucleotides. It cannot be admitted
as a coding-synonymous neighborhood by relabeling its constant reference
protein. Separate synonymous recoding validations in the paper are prior
evidence; they do not make the entire motif library a coding graph.

## Why intended sequences do not resolve processed products

The source already establishes that coding-region donors can remove the
reference stop and change the product. It also reports normal and
palindromic-insertion plasmid copies under the same reporter identities,
and ambiguity about the plasmid origin of some spacer-acceptor products.
Those findings are established prior work, not new discoveries here.
The intended FASTA alone cannot assign each recovered RNA to a full
realized plasmid genotype and translated product.

The source processing notebook filters to full-length/spliced categories,
applies barcode-count conditions and computes medians after concatenating
biological preparations. Its processed fractions must not be treated as
separate genotype-level biological replicates. Raw per-preparation files
exist, but a full sequence/junction-to-clone map and an explicitly chosen
new inference are needed before a processing-aware graph comparison.
The source code's treatment of multiple junctions requires further audit
if full transcript paths, rather than its reported category endpoint,
become the target. No source bug or corrective outcome is asserted here.

## Admission and stopping decision

Package C is an enabling source/identity result, not a new biological
headline. Recomputing published cryptic-splicing fractions would repeat
the source's main finding. No phenotype fit or graph reclassification is
launched until a distinct inference and authenticated realized mapping
exist. For package D, shared cell-line reporter labels could support a
regulatory-context question, but do not validate a codon-additive
synonymous model. For F, this source fails the changing-coding-genotype
criterion for its large motif library.

The concrete missing record is a sequence-linked full plasmid and
observed transcript/junction assignment for each assayed clone, with
barcode multiplicity and preparation identity preserved. Alternatively,
the new cPRIME coding source may supply a narrower processing-aware
application after its own genotype and amplicon gate. No outcome from
this Dao deposit has been promoted to fresh validation.

Primary evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC12297286/
and https://github.com/MustoeLab/Publications/tree/7a8779a484af797953a154966e8b06b4afced3ed/Dao_MPRAsplicing .
