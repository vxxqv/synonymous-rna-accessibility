# cPRIME local-state design and estimability audit

The selected design has only one additive constraint across 66 local states. Four of five exons are completely saturated by a model containing exon intercepts and recorded unary allele descriptors. The remaining constraint uses four states in exon 17. It is a balanced contrast between two intended-edit choices and two marker choices, not a conventional two-nucleotide-site factorial. This substantially limits any new interaction analysis.

No expression values, new source outcomes, phenotype missingness or model fits were inspected. The audit uses the previously fixed metadata selection. The separate proposed exploratory order-model protocol is ready for root review, but is not authorized for execution.

## Authenticated object and unknown coordinates

All 66 canonical IDs match a unique source-design target string. For every state, the two bracketed substitutions match the reported genomic reference bases, alternative bases and coordinate spacing. Every counted DNA word occurs exactly once in the alternative design context and includes both specified changes. Source reference strings agree wherever their designed windows overlap. These exact checks support local reference-aligned nucleotide identities.

They do not authenticate a complete realized haplotype or translated protein. The DNA words have different windows. `inventory_v1/coordinate_inventory.json` explicitly stores a nucleotide or null at every selected exon coordinate for each state. Sixty-two states have at least one unknown selected coordinate outside their own counted word. All other positions outside that word also remain unknown. No unknown was replaced with a wild-type base.

The model audited below is therefore a **recorded local-allele descriptor model**. Its zero indicator means that an allele is not one of the state's two recorded descriptors. It does not assert that a distant, unobserved coordinate is reference. The complete 66-state map and all exact reference alignments are preserved in `exact_v1/mapping.json`.

## Exact design results

| Exon | Local states | Intended-edit groups | Unary descriptor rank | Residual dimensions |
|---|---:|---:|---:|---:|
| 17 | 26 | 12 | 25 | 1 |
| 25 | 4 | 2 | 4 | 0 |
| 27 | 10 | 4 | 10 | 0 |
| 39 | 9 | 4 | 9 | 0 |
| 42 | 17 | 8 | 17 | 0 |
| Total | 66 | 30 | 65 | 1 |

There are 87 distinct recorded physical alleles and five exon-intercept columns. The resulting 66 x 92 matrix X has exact rank 65, a one-dimensional left nullspace and a 27-dimensional parameter nullspace. Fraction-based row reduction, explicit nullspace dot products, a modular lower bound of 65, and an independent unsigned-incidence graph formula agree.

For the latter check, treat physical alleles as vertices and recorded pairs as edges. There are 22 connected components, all bipartite. For one connected bipartite component, a vector annihilated by the unsigned incidence matrix must alternate sign along edges and hence has exactly one free value. Its rank is therefore number of vertices minus one. Summing gives 87-22=65. Exon intercepts add no rank because each is one half of the sum of its exon-specific allele columns. This is elementary known incidence-matrix algebra, not a new theorem.

Complementing each descriptor indicator within its exon and retaining the exon intercept gives the same exact column space and rank. This demonstrates coding invariance of the descriptor model. It does not convert an unobserved base into a measured reference allele.

The two-context model ranks are:

| Model of the latent descriptor score | Rows | Rank |
|---|---:|---:|
| Exactly shared coefficients across cells | 132 | 65 |
| Shared slopes, cell-specific exon offsets | 132 | 70 |
| Context-specific unary coefficients | 132 | 130 |
| Context-specific coefficients plus all observed-pair products | 132 | 132 |

The offset distinction is irrelevant to a model with a separate unrestricted nondecreasing response function in each exon/cell: its response function absorbs an additive offset. The latter three rank statements alone do not establish any biological fit.

## The sole additive constraint

Define four local states, preserving their exact genomic IDs:

| Symbol | Intended allele | Marker allele |
|---|---|---|
| A | chr17-31225100-A-T | chr17-31225106-A-G |
| B | chr17-31225100-A-T | chr17-31225109-T-C |
| C | chr17-31225112-C-T | chr17-31225106-A-G |
| D | chr17-31225112-C-T | chr17-31225109-T-C |

Every additive descriptor score obeys z_A-z_B-z_C+z_D=0, and this relation generates the entire left nullspace. The four counted words share an authenticated observed window from chr17:31225099 through 31225116 inclusive, covering all four changing positions. Their nucleotide Hamming-distance matrix in that common window is

```
0 2 2 4
2 0 4 2
2 4 0 2
4 2 2 0
```

Thus the two abstract factors are choices between different edit locations, not single changes at two nucleotide sites. No pair of these four states is a single-nucleotide neighbor. An independent exhaustive scan of all 17,667 within-exon four-state subsets found 129 whose common counted window covers every specified change, and zero complete two-site nucleotide factorials. This result is restricted to the selected 66 states and the required common observed window; it is not a claim that the full author library contains no useful controls.

Adding all 66 observed-pair products supplies an identity matrix on the observed states. It adds only one dimension to the pooled model, and one per cell in the context-specific model. Individual pair coefficients are not identifiable jointly with the unary columns. A fitted saturated pair model must not be reported as mechanistic evidence.

## Which contrasts are estimable

The 30 groups provide a basis of 36 observed within-intended-group contrasts. These compare combined local states carrying different markers. They are not measurements of an intended edit relative to a matched marker-free or wild-type genotype. Under the additive descriptor assumption they can identify differences between marker coefficients; without that assumption they remain conditional combined-state differences.

Absolute unary coefficients are not identifiable. A concrete null direction adds one to the coefficients of chr17-31225100-A-T, chr17-31225112-C-A and chr17-31225112-C-T, and subtracts one from chr17-31225103-T-C, chr17-31225106-A-G and chr17-31225109-T-C. It changes no selected-state score. The full right-null basis and the exact zero dot products are archived.

Only A, B, C and D have a within-context leave-one-state-out prediction that is estimable from the remaining rows under unrestricted descriptor additivity. For example z_D=z_B+z_C-z_A. The other 62 rows occur in no additive left-null relation, so their missing score is not fixed by the other 65 observed scores. No observed fit was used to reach this conclusion. A new prediction claim based on this same parameterization would therefore be largely underdetermined, not a credible broad validation exercise.

Two minimal counterexamples explain the limits. First, observing 00,10,01 with values 0,1,1 does not distinguish the additive missing value 11=2 from the interacting value 11=7. Second, latent rectangle values 0,1,2,3 satisfy the additive relation, but the strictly increasing transformation on these levels to 0,1,4,9 yields an observed-scale contrast of 4. A nonzero raw contrast therefore does not by itself exclude an additive latent score with a nonlinear response.

## Concrete next decision

Close broad descriptor-interaction fitting and full-genotype connectivity claims for this selected design. They would be saturated or unsupported by the measured sequence windows.

One bounded exploratory compatibility test is structurally capable: compare the actual shared-additive, context-specific additive and saturated observed-pair latent models using all selected states, both cells and the four previously declared interval-widening factors. `nested_protocol.md` defines the models, strict order constraints, exact rational witnesses and dual certificates. Only exon 17 can distinguish context-specific additivity from the saturated reference. The existing unrestricted shared-order comparison is weaker than shared additivity and must remain separately labeled. Root must review and freeze the new protocol before fitting already exposed B1 intervals.

For a conventional nucleotide interaction, the minimal missing evidence is an authenticated common-window WT, each single change and their double change, measured in the same assay with recoverable preparation identity and normalization. For a two-sided intended-edit/marker decomposition, the missing marker-free and intended-edit-free controls must be explicit. This audit does not propose new engineered sequences or assume those controls exist elsewhere.

## Source and prior boundary

The author sequence source is [Boettcher-Lab/cPRIME, commit 30125415304664c33370cd7801f820a345834e40](https://github.com/Boettcher-Lab/cPRIME/tree/30125415304664c33370cd7801f820a345834e40). Metadata were extracted from the local pinned source design, search-word inventory and Supplementary Data 8 identities. The current source paper is [El Kassem et al., Nature Communications, DOI 10.1038/s41467-026-78344-9](https://doi.org/10.1038/s41467-026-78344-9). Its published marker-removal and marker-substitution examples already establish that certain marker combinations affect RNA processing. Those published examples are exposed prior evidence; this audit does not claim their finding as new. The new Nature URL returned a cookie redirect during this bounded check; the already pinned local source text was available.

Linear estimability by row-space membership, equivalently orthogonality to a nullspace, is established machinery, described in the [official estimability documentation](https://search.r-project.org/CRAN/refmans/estimability/html/estimability-package.html). The standard four-state difference-in-differences and its reference dependence are also explicit in [Zhou and McCandlish, 2020](https://www.nature.com/articles/s41467-020-15512-5), equations 1 and 2. This audit contributes a source-specific design diagnosis and exact finite certificates, not new generic estimability or epistasis theory.

## Reproduction and inspection limits

All three metadata audit scripts and the synthetic ordinal preflight use the Python standard library. Run from the workspace root, choosing new output directories:

```
python work/research_cycle8/context_design_audit/audit_design.py --root . --out NEW_EXACT_DIRECTORY
python work/research_cycle8/context_design_audit/independent_checks.py --inputs NEW_EXACT_DIRECTORY --out NEW_CHECK_DIRECTORY
python work/research_cycle8/context_design_audit/compile_design_inventory.py --inputs NEW_EXACT_DIRECTORY --out NEW_INVENTORY_DIRECTORY
python work/research_cycle8/context_design_audit/ordinal_preflight.py --out NEW_SYNTHETIC_DIRECTORY
```

The completed run verified 1,874 exact nullspace dot products, the modular rank sandwich, all quartet comparisons, reference recoding, and eight invented ordinal-contract cases. A second algorithm and representation were used by the same agent with inherited project context; no blinded, external or separately authored review is claimed. The complete hashes, runtime, input mapping and file inventory are retained. No raw workbook, count table, new biological result, or phenotype model was accessed by these scripts.
