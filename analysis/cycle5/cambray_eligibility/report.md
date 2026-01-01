# Cambray joint-function eligibility decision

## Decision

**Conditional go for a small, explicitly descriptive joint-trait accessibility study. No phenotype analysis is authorized by this report.** The deposited design supports 47 complete synonymous binary squares, comprising 181 designs, 42 peptides, 43 observed synonymous components and 32 design seeds. All squares remain after excluding the entire previously exposed seed-1 family. Whether their required phenotype entries are present and usable remains unexamined. A complete unrestricted synonymous landscape, independently paired assay preparations, calibrated joint uncertainty, and physical full-backbone identity are not established.

The most discriminating next question is whether two acceptable endpoints have a protein-preserving route and a growth-preserving route, but no route preserving both on the same declared graph. This is different from a correlation, an expression-growth tradeoff, or a Pareto frontier. It is also different from evidence for epistasis: opposing, exactly additive traits can produce the proposed obstruction. A positive observation would establish a restriction imposed by the chosen pair of assay requirements, subject to the stated graph and measurement scope.

## Scope and object filter

This audit covers the authenticated Cambray author deposit, its original article and supplementary methods, closely related primary studies, and all deposited design identities. Source figures and published aggregate results are exposed literature, not untouched validation outcomes. The code reads only the first five identity fields from each integrated-data row. New phenotype values are neither decoded nor evaluated. The extracted sequence inventory contains no outcome columns. General codon-use reviews, unrelated clinical or viral papers and broad predictive feature studies were excluded from the decision evidence. The targeted backward and forward search is not an exhaustive novelty review or evidence that no relevant prior exists.

## Exact sequence evidence

`graph_summary.json`, `synonymous_edges.tsv`, `binary_faces.jsonl`, `annotated_squares.json` and `square_dependencies.json` give the inspectable inventories. The full ZIP is 130,122,456 bytes with SHA256 `2758ed63e469c77b42c01e3c022452ef52ebf8d1bd4a86bd4c5b502e44b7bbc0`. Its initial bytes exactly match the earlier cached prefix. The sequence-only TSV SHA256 is `9fcd79a97e3478c68411dee6cb853ad77e249d1d93b11bc4af58d2fcc44d4e7e`.

| Quantity | Exact result | Interpretation |
|---|---:|---|
| Unique design IDs and 96-nt sequences | 244,000 each | No exact sequence duplicates |
| Intended 32-aa peptides | 157,278 | This is not one synonymous library |
| Internal stop codons | 0 | Standard code, stated frame |
| Design seeds | 56 | Fifty-five have 4,374 designs; seed136 has 3,430 |
| Synonymous single-nucleotide edges | 20,031 | All remain within their seed |
| Vertices incident to synonymous edges | 34,426 | Most designs are isolated in the sampled graph |
| Two-step routes with Hamming-two endpoints | 6,085 | 5,991 distinct endpoint pairs; phenotype coverage unknown |
| Largest observed synonymous component | 23 vertices | No claim about unmeasured vertices |
| Complete binary squares | 47 | Two distinct codons vary in every square |
| Complete binary cubes | 0 | No complete three-coordinate factorial in this inventory |
| Squares that equal an entire observed component | 10 | The other37 require checking observed detours |
| Squares after conservative prior-exposure exclusion | 47 | All181 constituent designs survive |

The label `rep` identifies three nearby designed sequences, not experimental replication. The missing944 third designs in seed136 are explicitly documented in the original design methods, not a newly inferred assay dropout. The apparent completeness of 1,749 peptide groups as products of their *observed* alleles does not make them complete products of all synonymous codons; 1,747 of these products have only two vertices.

The original codon table and an independent standard-code table agree at all64 codons. A separate algorithm then compares all3,916,771 pairs within translated-peptide groups as strings, reconstructs cycles instead of expanding corners, and exactly reproduces the full edge set, all47 squares, absence of cubes and component sizes. This is algorithmic corroboration by the same analyst, not a separate-person review. It does not independently re-enumerate nonsynonymous edges. See `graph_corroboration.json`.

## Coding frame and full context

The original methods place the variable96nt directly after a fixed ATG and before a common linker, TEV site and GFP reporter. Supplementary Note1 explicitly states that the start codon is fixed and the designed segment encodes32 codons. Author code supplies the upstream sequence ending in ATG. The original oligonucleotide design gives the fixed flanks and KpnI/BamHI insertion scheme. Addgene113408/113409 identify the shared pGC4750 backbone and C-terminal linker/TEV/GFP context.

Consequently, equality of the translated96nt segment proves equality of the **intended full reporter ORF** under the authenticated common, codon-aligned context: for any common in-frame prefixP and suffixS, translation of P+u+S equals that of P+v+S when u,v are96nt and have identical32-aa translations. The exact common suffix is not needed for this implication. This does not prove that every physical clone had an otherwise identical genome or plasmid; source methods acknowledge mutations outside the sequenced region and remove conspicuously low protein observations on that basis. Full Addgene sequence files require login. No restricted download was attempted. The public design and common-context evidence support intended synonymy, with this physical-genotype limitation explicit.

## Measurement contract

The deposited rows align traits by designID. They do not align protein and growth measurements on the same culture or produce paired cross-assay biological replicates. The following semantics are based on the original methods, Supplementary Note1 and the saved author processing scripts, not on new observed values.

| Endpoint | Preparations and units | Processing and dependence | Permitted interpretation |
|---|---|---|---|
| Normal-initiation protein, `clean.lin.prot.*` | Four independent biological FACS preparations; two instruments; transformed fluorescence in arbitrary units scaled1 to100 | Modal-bin selection, replicate rescaling, low-outlier removal, quantile normalization and a shared nonlinear scale; counts and bins are not additional biological replicates | Processed reporter output in rich MOPS, not protein molecules, benefit or fitness |
| Normal-initiation growth, `growth.R60one/two/three` | Three parallel competition lines derived from one initial library culture; relative growth over approximately60 generations | Shared initial abundance reference, library-compositional normalization and generation calibration; low counts removed | Relative pooled competition phenotype in that environment, not absolute doubling time |
| `growth.aggregate.R` | Combination of the R60 mean with R13 and R28 from a fourth follow-up experiment | Geometric aggregation across times and preparations; R13 and R28 are not independent biological replicates | Author integrated growth point estimate; retain separate preparation/time sensitivity |
| RNA abundance/decay | Two independently cultured samples, transcription arrest, six times0 to40min | Shared mean DNA denominator, spike-in correction, empirical time shift for preparation2, nonlinear protected-component fit, quality and sign filtering | RNA abundance ratio and fitted decay summaries; neither transcript count nor directly calibrated uncertainty |
| Polysome load | One recovered profile split into fractions, not five biological replicates | Fractions6/7 discarded; firstfive normalized within construct; weighted mean requires at leastfour finite fractions | Relative polysome distribution, not independent replicated translation rate |

The strongest initial pairing is normal-initiation protein plus growth in normal rich MOPS, because its condition labels agree. Facilitated-initiation protein must be paired with facilitated growth, and minimal-medium growth is a different environment. RNA lifetime need not be larger to be better. GFP output likewise is an assay objective and can impose a growth burden. Any eventual bounds must be justified as reporter-output requirements or clearly labeled operational ranges; arbitrary percentile bounds cannot establish natural biological acceptability.

The author scripts make calibrated replay a separate gate. Protein prose and Note1 specify four consecutive fluorescence bins, while the regular-initiation code calls a three-bin filter with a different depth setting. RNA code assigns the first replicate twice and later uses an undefined second object, so the deposited fragment cannot establish an executable generating pipeline by itself. Growth code uses aliases and intermediate objects absent from the fragment. These are provenance/replay discrepancies, not evidence that the published measurements are false. The deposited processed values can support a clearly identified descriptive analysis without pretending that a repaired script reproduces the author's uncertainty.

## Identifiability before outcomes

For square endpoints00 and11, suppose both endpoints pass both requirements. Let a,b be pass indicators for the two intermediate vertices under trait1, and c,d the corresponding indicators under trait2. Joint-only obstruction is exactly

`(a or b) and (c or d) and not ((a and c) or (b and d))`.

It holds only for the two opposite-intermediate patterns. `square_design_check.py` checks all16 patterns against graph search. Its positive control uses trait1 values `[1,0,2,1]` and trait2 values `[1,2,0,1]`, threshold1, with both endpoints retained. Each field is exactly additive, each trait has a path, and their intersection has no path. Aligned or identical acceptability sets, an absent individually permissive path, or a jointly unacceptable endpoint are negative controls. A complete square has one residual contrast beyond its three-parameter additive model, unlike the prior saturated three-state design, but observing a joint barrier alone does not test that contrast.

The square obstruction concerns paths confined to its two changing positions. Before making any broader observed-component claim, test every deposited same-peptide detour in its component. A path found in this sampled graph is a witness. A missing path proves only sampled-graph disconnection; unmeasured synonymous routes remain unknown. Forty-seven overlapping squares are neither47 independent biological replicates nor47 independent design families.

For an eventual frozen analysis, retain all94 opposite-corner pairs and all failures; specify missingness before evaluating signs. Evaluate an explicitly defined threshold domain rather than select favorable thresholds. Retain the entire prepared vector when assessing preparation sensitivity, so common normalizations and shared mutants remain coupled. Cross-products of four protein and three growth preparations are dependent sensitivity scenarios, not12 independently paired experiments. Threshold or preparation ranges are not confidence intervals. Calibrated joint confidence claims remain gated on a defensible measurement model or raw-data replay.

## Stop and continuation conditions

Proceed only to a root-approved protocol and selective finite-coverage stage after this literature/eligibility decision. Freeze the sequence-selected units, both diagonals, boundary convention, condition pairing, admissible bounds, missingness treatment, primary graph scope and all negative controls before new phenotype values are read. The exact prior503IDs are listed in `previously_exposed_ids.json`; exclude the whole seed1 family and its peptide identities from any nominal untouched evaluation. Published figures and source aggregate exposures also remain disclosed.

Stop the biological claim if there are no complete paired-trait units, no interpretable requirement ranges, no stable obstruction across declared sensitivity scenarios, or no evidence that separates the proposed graph statement from threshold choice and sampling. A null result in a well-covered declared domain is informative. If only processed means are usable, the strongest claim is descriptive accessibility in an observed reporter-design graph. It is not replicated mechanism, natural selection, a complete sequence-space barrier, or experimental improvement from engineered interventions.

## Deliverables and reproducibility

Run `recover_sequences.py` only into a fresh audit location to reproduce selective extraction; the author ZIP is already retained. `audit_graph.py` requires fresh graph output paths. `check_graph.py` is a distinct exact corroboration. `square_design_check.py` produces family counts and the synthetic truth-table controls. The final provenance manifest hashes scripts, source records and result inventories. No frozen project result, manuscript or release was changed.
