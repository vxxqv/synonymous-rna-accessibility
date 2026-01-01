# Closest-prior and discrimination map

Checked 2026-10-25. This is a bounded original-source review, not a global novelty certificate. The proposed proof is a formal coding/product corollary of classical vector rearrangement. The biological usefulness would be a falsifiable compromise bound under explicit additive/local-intermediate assumptions, not a new bottleneck algorithm or evidence that any measured trait is additive.

## 1. Constant trait count is not enough for exact tractability

**Jean Cardinal, Erik D. Demaine, David Eppstein, Robert A. Hearn and Andrew Winslow.** *Reconfiguration of satisfying assignments and subset sums: Easy to find, hard to connect*. Journal version: Theoretical Computer Science806(2020),332-343, DOI [10.1016/j.tcs.2019.05.028](https://doi.org/10.1016/j.tcs.2019.05.028). Exact source read: [arXiv1805.04055v1, section4, Theorem4.7 and proof](https://arxiv.org/html/1805.04055v1); the [author-hosted conference PDF](https://andrewwinslow.com/papers/reconfig-cocoon18.pdf) has the shorter statement as Theorem4.

The source proves PSPACE-complete hypercube reachability constrained by a constant number of linear inequalities. Its displayed full reduction uses70 dimensions: two value dimensions and one count dimension per each of23colors, plus one key dimension. The relevant one-sided requirements are46lower bounds and24upper bounds, not a two-trait hardness theorem. Coefficients and ranges are polynomially bounded in the reduction, so their bit lengths are logarithmic. This is load-bearing negative prior: the proposed R=0 decision on repeated twofold Lys codons is the same problem after translating variables into AAA/AAG wobble choices. Its trait rank is at most70. That translation is straightforward, not a new generic complexity theorem.

Overlap label: **known after translation of notation**. Do not infer tractability from fixed trait count, fixed effect rank or absence of epistasis.

## 2. One scalar interval already has a relevant hard boundary

**Takehiro Ito and Erik D. Demaine.** *Approximability of the subset sum reconfiguration problem*. J.Combin.Optim.28(2014),639-654, DOI [10.1007/s10878-012-9562-z](https://doi.org/10.1007/s10878-012-9562-z). [Author paper and current correction notice](https://erikdemaine.org/papers/SubsetSumReconfiguration_JOCO/), [original manuscript](https://erikdemaine.org/papers/SubsetSumReconfiguration_JOCO/paper.pdf), [two-page erratum](https://erikdemaine.org/papers/SubsetSumReconfiguration_JOCO/erratum.pdf) read.

The configuration is an item subset with its scalar sum between a lower threshold and an upper capacity; one item is added or removed per step. This is one interval trait, or two oppositely oriented lower-bound traits, on a cube. The correction replaces Theorem1 and Corollary1 with NP-hardness, withdrawing their stronger NP-hardness claim. The original also studies maxmin compromise, approximation, and conflict graphs. Its conflict-graph PSPACE result is adjacent to the existing pairing/independent-set representation, but not needed for C5.1. Do not cite the original abstract's uncorrected strong-NP wording. The new bound is not a resolution of the exact two-trait PSPACE question discussed in Cardinal's source.

Overlap label: **known after translation of notation** for the interval reconfiguration formulation; **partial or adjacent result only** for our rank/local-codon bound.

## 3. The vector-ordering engine is classical, with current refinements

**V.S.Grinberg and S.V.Sevastyanov**, *Value of the Steinitz constant*, Funct.Anal.Appl.14(1980),125-126, DOI [10.1007/BF01086559](https://doi.org/10.1007/BF01086559). Original metadata authenticated through [MathNet](https://www.mathnet.ru/eng/faa1805); full original text was not reliably accessible, so no original-language theorem translation is claimed.

The exact external statement used is verified in **Imre Barany**, *A matrix version of the Steinitz lemma*, J.reine angew.Math.809(2024),261-267, DOI [10.1515/crelle-2024-0008](https://doi.org/10.1515/crelle-2024-0008), [author manuscript arXiv2308.10102v2, Theorem1.2, page2](https://real.mtak.hu/212570/1/2308.10102v2.pdf). Zero-sum vectors in a norm unit ball in dimensiond can be reordered with partial sums bounded byd; the paper attributes that improvement to Grinberg-Sevastyanov. Restricting the ambient infinity norm to the endpoint-effect span gives r, not raw trait count. Our collinear-dummy argument and factor-path lift are explicit applications of that result. The local endpoint-box excursion is the part required to respect actual codon mutation edges.

Overlap label: **formal corollary not stated** for C5.1, conservatively scoped to the checked sources. No priority claim for rank reduction, dummy augmentation or rearrangement.

**Gergely Ambrus and Rainie Heck**, *A note on the Steinitz lemma*, Mathematika72(2026),e70085, DOI [10.1112/mtk.70085](https://doi.org/10.1112/mtk.70085), publisher full text read, particularly framework/definitions and historical discussion. It studies non-zero-sum variants and their relationship with classical Steinitz constants. This makes a novelty claim for the generic endpoint-segment ordering particularly inappropriate. Our conservative r constant is not advertised as the sharpest known bound. New matrix/colorful Steinitz refinements are not imported without checking their different ordering contracts.

## 4. Joint-function pathways already overlap strongly with Kirsch2026

**Orson Kirsch, Nicole Wood, Steven A.Redford and Kabir Husain.** *The fitness landscape of overlapping genes*. [arXiv2604.00602v1](https://arxiv.org/html/2604.00602v1), submitted1April2026; explicitly a preprint. Original ResultsII.1-II.2 and MethodsIV.6 read.

It examines exact raw codon compatibility and then joint functional preferences inferred using two protein-family Potts models. Amino-acid sequences can change while retaining inferred family compatibility. Its mutational-path search orders the differing endpoint nucleotides, using a genetic algorithm and a maximum along-path sum of absolute energy z-scores; internal stops are forbidden. Thus joint functions, frustration, Pareto tradeoffs and near-neutral paths are major conceptual overlap, not an unoccupied topic. Our graph fixes specified peptide(s), our phenotype constraints can use separately justified assay intervals, and C5.1 assumes a product/additive approximation. The formal bound and explicit codon-intermediate correction differ from that heuristic Potts study; neither is established as a new general principle. Its permutation search does not certify absence of a route using temporary background changes, as our Q3 fixture illustrates. We do not reinterpret inferred family scores as measured trait preservation.

Overlap label: **partial or adjacent result only** for C5.1; substantial conceptual overlap for the proposed biological question.

## 5. Generic path, convexity and spectral machinery is not the contribution

**Mehryar Mohri**, *Semiring frameworks and algorithms for shortest-distance problems*, J.Autom.Lang.Comb.7(2002),321-350, DOI [10.25596/jalc-2002-321](https://doi.org/10.25596/jalc-2002-321), [author manuscript](https://cs.nyu.edu/~mohri/pub/jalc.pdf) read for its algebraic shortest-distance framework. Min/max path recurrences and a vertex-violation scalarization are standard. The present R is a standard bottleneck quantity. This is not a new multiobjective-path algorithm; the biological specification of scales and acceptable ranges supplies the interpretation, not generic algorithmic novelty.

Graph-convex acceptable sets have convex intersections by definition. An additive superlevel on a cube need not be graph-convex: for x1+x2>=1, the shortest route10-00-01 is forbidden, although another shortest route exists. Full codon factors can have triangles, so invoking a binary median-graph theorem without a verified quotient/lift is invalid. Existing project Gopalan/Schwerdtfeger and multiple-target RNA design comparisons remain relevant, but none is a new external leaf needed for this proof.

For a filtered graph A, its combinatorial Laplacian has nullity equal to its component count: z^T L_A z=sum_{uv in E_A}(z_u-z_v)^2. This standard identity is proved directly and has no novelty. A principal submatrix of the base Laplacian is generally a grounded operator, not this induced Laplacian: Q2 restricted to00,11 gives principal diag(2,2) but induced zero2x2. Base/product spectra therefore do not decide joint acceptability. Effective resistance can detect disconnection only after the correct filtered graph is known. Eigenvalues may describe robustness or mixing of a measured graph, but would not by themselves identify the restrictive biological function. A spectral formula alone is rejected as a contribution.

The [official Fan Chung spectral graph lecture archive](https://cseweb.ucsd.edu/~fan/research/spectral_graph.html) and its first-chapter/full-book links were located, but full-text retrieval timed out. No external theorem verification is claimed from that inaccessible text. The identity above is self-contained.

## 6. Relation to the project's existing uncertainty results

For fixed finitely many trait constraints, requiring every trait at every visited vertex is the fixed-route robust quantifier already treated in the project: set one scenario score to minus each scaled violation. Then -R=max_P min_k min_{v in P}(-psi_k(v)). The biological labels differ from uncertainty scenarios, but the algebra is the same. No new quantifier theorem is claimed.

The genuinely useful discrimination is narrower: the product/additive model guarantees a finite compromise bound, with a mandatory codon-intermediate correction and a deterministic2epsilon approximation allowance. An observed violation of a calibrated bound would reject at least one model/measurement assumption, rather than automatically proving a folding mechanism. Uniform error over every relevant genotype is not supplied by in-sample explained variance, mean prediction error or a small omitted singular value. A low-rank approximation only enters via C5.3 if its entire-state error is bounded externally or deterministically from complete contribution tables.

## Search and stopping record

Queries covered exact and synonymous combinations of joint/multiobjective paths, constrained hypercube paths, knapsack/subset-sum reconfiguration, graph-convex/median preferences, Steinitz arbitrary/nonzero sums, vector sequencing/scheduling, and Kirsch's exact title. Searches used the web index and primary publisher/arXiv/ECCC/author-hosted sources; backward references led from Cardinal to Ito-Demaine and its erratum, and from2026Steinitz work to1980/2024forms. The second strategy checked coding and synonymous terminology against vector-rearrangement vocabulary. No close coding-specific Steinitz formulation was returned, but that absence is not a novelty result. Relevant Russian scheduling/strip literature was located at metadata level, not fully checked. No claim of exhaustive coverage, first use, sharp general constant, or global priority is justified.

Recommendation: send C5.1-C5.3, their counterexamples and the exact finite contract to an independent reviewer. Integrate only if the bound materially improves the empirical question and its assumptions can be made explicit. Otherwise retain it as a negative prior/discrimination map, rather than expanding generic mathematics merely to enlarge the paper.
