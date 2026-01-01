# Independent FAS sequence and graph verification

**Verdict: PASS for the exact declared design and native included-isoform coding claims. Measured completeness, actual minigene translation and protein output remain unverified.** The reconstructed ancestor used as `WTSEQ` by the author pipeline is not in the native human protein class.

This review was separately authored with inherited project context. The other lane's source and reported counts were read, so it is not a blinded check. Its enumeration was not executed or imported. The independent implementation uses a literal standard-code table, a restricted explicit GenBank parser, codon-level product enumeration, all-pairs Hamming adjacency, and breadth-first graph distances. It does not use Biopython or an author numerical module. The finite calculation took about0.4s and used no phenotype table or count file.

## Exact contract and source identity

The question is about the author's existing declared63nt design under the native complete included-isoform CDS, not a newly designed sequence set. The native record is NM_000043.6, protein NP_000034.1, SHA256 `e8964c4db7d19acd5f6a1075106ec67087ad4bc5a8dc3a986f7eccb0a7f31ca7`. Independent translation of the entire CDS agrees with its annotated335-residue product and terminal stop.

The cached author filter has exactly the same design regex for the input and output loops. Its content has Git blob SHA1 `ba3858d75e4072c723c53c4842a9f88a792aa490`. Independent reconstruction of all52 cached Git tree objects links that blob to canonical tree `b526e2a0fa5c1b72be91b8176677a83d9b57bf28` in author commit `b87064c39a553b0759d86c2b428813c58dc17ae2`. The inventory's top-level sha denotes the requested commit, not the canonical tree; the explicit object reconstruction avoids treating those as equal. This authenticates the cached source chain without a new network retrieval. It does not establish that every declared sequence survived the experimental processing filter.

## Coding phase and complete enumeration

The63nt design matches exactly one annotated exon, transcript positions585..647 inclusive. The native CDS is80..1087 inclusive. Thus the exon begins at zero-based CDS offset505, or codon position2. Its first codon needs one preceding CDS nucleotide and its last codon needs two following CDS nucleotides. Translating the isolated exon from its first base would use the wrong phase. The verifier retains both flanks and substitutes only the existing exon.

There are11 variable nucleotide positions: one three-state and ten binary domains. This gives `3*2^10=3072` distinct designs. Local codon amino-acid alternatives give80 distinct complete native-context protein products. An independent product of local amino-acid multiplicities agrees with direct full-CDS enumeration. The protein-class size histogram is:

| Class size | Number of classes |
|---:|---:|
|8|16|
|16|4|
|24|32|
|48|8|
|72|16|
|144|4|

The weighted sum is3072, no design is duplicated, and no enumerated protein class contains a premature stop. All72 native-reference class strings agree exactly with the previously reported sequence list. This is a finite design check, not an observation of72 measured constructs.

## Why the human reference class is a Cartesian product

Equal native protein requires each affected codon independently to retain the native amino acid, because the CDS is single-frame and all other nucleotides are fixed. The surviving alternatives at the five nontrivial codons are:

| Native CDS codon number | Existing allowed synonymous codons | Local nucleotide graph |
|---:|---|---|
|177, leucine|CTT, CTG, TTG|CTT - CTG - TTG, a three-vertex path|
|179, leucine|CTC, CTT|one binary edge|
|180, leucine|CTG, CTT|one binary edge|
|182, leucine|CTC, CTG, TTG|CTC - CTG - TTG, a three-vertex path|
|183, proline|CCA, CCG|one binary edge|

Other affected codons have one surviving native-class alternative. The full native class is therefore exactly `P3 x P3 x K2 x K2 x K2`, or `P3 x P3 x Q3`. This follows from disjoint codon coordinates, not from fitting a graph name to numerical totals. The direct72-string graph agrees with the Cartesian adjacency criterion on all2556 unordered pairs.

The vertex count is `3*3*2^3=72`. The edge count is `2*(2*3*8)+3*(3*3*4)=204`. Its degree histogram is32 vertices of degree5,32 of degree6 and8 of degree7. Independent breadth-first search confirms connectivity, diameter7 and graph distance equal to nucleotide Hamming distance for all5184 ordered endpoint pairs, including identical endpoints. Native-reference distance shells have sizes1,5,12,18,18,12,5,1.

## Material qualification about the reference

The author count script calls a reconstructed ancestor `WTSEQ`. Substituting that exon into the same native coding context does **not** preserve the native human protein. Therefore it is false to describe all72 human-reference variants as synonymous relative to the author's ancestor. The correct object is the class preserving the human native included-isoform protein. A later assay normalization against the ancestor, if used, must be described as such and cannot supply protein identity. The full3072-design library also is not one synonymous library: it spans80 native-context proteins.

Even the72-member statement concerns the native included isoform. The experiment is a minigene RNA splicing assay. This verification does not establish that its transcripts are translated in the same frame, that a common native protein is produced experimentally, or that changing exon inclusion preserves all protein products. It supports the phrase “variants preserving the native included-isoform coding sequence product,” not “measured protein-preserving biological routes.” Exact construct context and assay metadata remain separate gates.

## Exposure incident and limits

A source read of `source_author_source_0.txt` displayed its five-row illustrative genotype/Occurrences table. Those are source README examples, not a count file, but they are conservatively recorded as exposure. No occurrence value was used in enumeration, classification or graph construction. `provenance_and_exposure.json` records the five exact sequence IDs and source lines without reextracting their numeric values; exactly one is in the72-member human reference class, the native reference itself. Root was notified immediately. Do not describe this review as wholly outcome-unexposed or silently use that native example as untouched coverage evidence. No splicing score, experimental count table or outcome-dependent missingness was inspected.

The design filter is an unanchored shell grep, so declared regex membership is not by itself authentication of experimental sequence length or read parsing. This audit matches the63nt pattern in full when enumerating and locates it uniquely in the native annotation. Future observed-data readers must separately enforce exact length/alphabet/design identity and authenticate the actual count schema; none was run here. Pinned current RefSeq identity is not a demonstration that every physical minigene clone was exactly the intended sequence.

## Reproduction and receipt

Run `work/env/Scripts/python.exe work/research_cycle6/fas_independent_verification/verify.py` in a fresh copy without existing result files. `check_provenance.py` verifies the cached Git objects and records the source-example exposure. Both refuse to overwrite their outputs. No author code is executed.

`manifest.json` records the exact inputs before and after, interpreter, command, current repository revision and dirty state, finite domain and result hash. The Mathbox v2 manifest validator passes against the workspace. A first manifest serialized declared result paths as objects rather than strings; its script, unchanged numerical result and invalid manifest are retained in `history_manifest_schema/`. Only that serialization field changed before the final rerun, whose numerical result is byte-identical. This was a receipt-format correction, not a sequence or scientific change.

The computation result SHA256 is `12c4f87cc30a8075404be708659837cc75df3f28908af55042f925da1841141f`. No graph or coding counterexample was found within this complete declared domain. The source-reference distinction and experimental-translation limits are consequential and must accompany any biological follow-up.
