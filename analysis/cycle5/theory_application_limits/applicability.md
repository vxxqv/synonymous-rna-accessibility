# Application limits of the reviewed rank/codon bound

Verdict: **applicable as a conditional theoretical interpretation and a posteriori finite-square check, not as an independently predictive or mechanistic headline on this source.** The available graph metadata do not establish the additive/full-product/global-error hypotheses for the larger sampled components. No phenotype values were opened, and no empirical rank, effect size or residual was fitted.

The source C5 proof is the reviewed `joint_math/claim_and_proof.md`, SHA256 `d5ef26923060fd587bb35b0c2d03d6ae39c34fda1e9a967f87e723e00d5cadb1`. Its separately authored review is `independent_math_review/review.md`, SHA256 `559036a4f3571ce2aeadebfe8e59f49773b45b00aa7c770ca003c288acc90d37`. That review approves the product/additive/closed-box result subject to inherited hypotheses and explicit integration clarifications. It does not validate biological additivity or uniform error. This report is a focused applicability audit by the original proof author with inherited context, not a fresh independent re-proof or external referee endorsement.

## Applicability table

| Object or claim | What is established without opening outcomes | What would be needed for stronger use | Permitted role |
|---|---|---|---|
| One of the47 complete squares | Exactly four intended-synonymous genotypes, two distinct codons, two binary single-nucleotide factors; the declared Q2 is a complete product | All required measurements finite, and separate measurement-error assumptions for latent claims | Exact finite-square graph analysis; eta=0 on this declared Q2 |
| Four-point additive projection | A rank3 design and one residual contrast per trait; the observed sup residual is exactly absolute contrast/4 | No additional assumption for this observed four-point identity | Descriptive contrast and identity check, only after prediction sealing |
| Approximation bound on observed Q2 | Once all four scores and scales are known, their projection satisfies a uniform residual bound over these four observed points | None for the algebra; biological interpretation still requires declared meaningful scales/requirements | A posteriori bound on that observed field, not independent model validation |
| True latent scores on the same Q2 | The measured projection alone establishes no true-score error bound | An independently justified uniform bound between the true scaled traits and this additive model at all four points | Conditional robustness statement only |
| One withheld corner | Three points fit a three-parameter additive model exactly | A target-independent deterministic error certificate or genuine held-out empirical evaluation | Frozen three-known-corner benchmark; no error certificate from saturated training fit |
| Larger observed components | The43 selected components have289 vertices; ten squares equal their whole observed component, while37 have observed detours | Exact product decomposition of the declared component plus an additive representation or uniform error over every vertex of that product | Direct graph connectivity and missing-state bounds; no automatic C5 extension |
| Complete synonymous component | Standard-frame codon factors are a mathematical product after taking their actual connected components | Contributions at all factor states, local codon-path eta, and a globally valid additive/error model on that declared product | A conditional theoretical consequence, not supplied by the sampled source |
| Low exact empirical rank | At most two changed-factor columns on Q2, regardless of trait count | A justified model of intrinsic rank rather than treating noisy exact rank as latent rank | Model description with explicit assumptions; no biological redundancy inference |

## Exact four-point identity and what it does not estimate

Use corner order(00,01,10,11), h=(1,-1,-1,1), and the observed trait vector y. The additive design spans the constant and two binary coordinate columns, whose orthogonal complement is span(h). Since h dot h=4,

```text
C = y00-y01-y10+y11,
y_add = y-(C/4)h,
max_v |y(v)-y_add(v)| = |C|/4.
```

For positive fixed trait scales a_k, the scaled uniform residual is exactly

```text
epsilon_observed = max_k |C_k|/(4 a_k).
```

The projected two factor effects are

```text
v1 = (Y10+Y11-Y00-Y01)/2,
v2 = (Y01+Y11-Y00-Y10)/2,
```

where Y denotes the scaled vector field and v1 changes the first binary coordinate. Their rank is at most2 and M=max(||v1||_infinity,||v2||_infinity). On the declared binary product the factor moves are direct edges, so no unmeasured local intermediate is traversed and eta=0. Applying the reviewed approximation corollary to **this four-vertex graph and these observed values** gives R_observed<=rM+2epsilon_observed.

This is a valid deterministic conclusion, but choosing the approximation allowance from the complete observed residual makes every four-point field satisfy the approximation hypothesis by construction. It imposes no independently testable model restriction on that same field. A violation of the resulting mathematical bound would indicate a calculation or hypothesis-translation error, rather than constitute independent evidence for a restrictive biological mechanism. R_observed itself is already exactly

```text
R_observed = min(psi(01),psi(10))
```

when endpoints00 and11 are acceptable. The two candidate simple paths are explicit. On this domain the rank bound therefore adds neither an efficient search advantage nor independent evidence for its model. It can be loose, or can illustrate sharpness. The exact field(0,-2e,-2e,0) in one trait has a constant additive projection-e, rank0, M=0 and uniform observed residual e; its endpoint-relative R is2e, attaining the corollary. All four values were used to calculate the certificate, so this tight identity is still nonpredictive.

The four-point residual is not the measurement error relative to a latent phenotype. A perfectly flat observed field has zero contrast even if its true intermediates are arbitrarily low because of uncalibrated vertex-specific error. Nor does the residual control unmeasured genotypes. On a three-vertex path with only endpoint values observed as0, a constant0 additive model fits those observations exactly, while an unmeasured middle score-L yields R=L for arbitrarily large L. Average held-out error, a training residual or preparation scatter does not supply the missing uniform guarantee.

There is a useful domain distinction: an upper-bound witness confined to a prespecified measured Q2 need not assume errors are bounded on unvisited unmeasured codons. Apply the corollary to that declared Q2, then use that its routes also lie in the larger base graph. This can bound the larger graph's R from above conditional on a valid true-score error bound on the square. It cannot establish that a sampled disconnection persists after adding unmeasured routes, and it does not authenticate a full-graph additive model.

## Exact rank is fragile to measurement noise

Start with the exactly additive opposing traits F=(x1-x2,-x1+x2). Its effect columns(1,-1) and(-1,1) have rank1. Perturb only trait2 at vertex11 by an arbitrary nonzero rational d. The four-point projection then has effect columns

```text
v1=(1,-1+d/2), v2=(-1,1+d/2), det[v1 v2]=d.
```

Thus the exact rank is2 for every nonzero d, however small. Numeric full rank of measured effects is not evidence for two independent latent mechanisms, and rank1 is not a stable point estimate under arbitrary noise. Treating decimal or binary floating values as exact rationals only makes the arithmetic rank exact for those stored numbers. It does not make the biological rank exact. A truncated singular-value rank is a different modeling choice; any discarded component must be included in a separately justified uniform approximation allowance over the declared domain. A post hoc numerical tolerance is not such a calibration. Invertible positive trait rescaling preserves exact rank but changes M and numerical conditioning.

## Product coverage and missing codons

The sequence-only inventory contains no complete three-dimensional binary cube. That fact does not prevent small nonbinary products, but it supplies no large factorial domain on which to test the bound. An observed same-peptide component is an induced sampled graph, not automatically the product of its observed codon choices. Even a product of all observed choices need not include all synonymous codons. Serine component gates and Leu/Arg intermediary codons still matter on a complete standard-code domain; eta=0 for the two directly connected binary square factors must not be reused for those larger factors.

A damaging sampling example is the induced Q4 path

```text
0000,0010,0011,1011,1111
```

with the exactly additive scalar score f=x1+x2-x3-x4 and threshold0. Its scores are0,-1,-2,-1,0, so the sampled graph has R=2. Naively treating its four observed coordinate choices as independent factors gives rank1, M=1 and eta=0, hence the invalid sampled claim R<=1. The actual full Q4 admits the route0000,1000,1100,1110,1111 with no negative score, so R_full=0. This is missing product coverage, not epistasis or an error in the reviewed theorem. In general R_full<=R_sampled; an upper bound on R_full does not upper-bound R_sampled.

The valid local square witness transfers to a larger graph as an existence witness. A local square obstruction transfers in the opposite direction only after ruling out every additional route, which the missing sequence space prevents. The protocol correctly requires testing all observed component detours and preserving unknown-state qualifications.

## A posteriori projection versus held-out prediction

Any three corners specify a unique additive completion; their training residual is zero because that fit is saturated. The fourth value can be any real T without violating those three observations. Before unsealing, neither its contrast nor epsilon_observed can be calculated without target leakage. M and rank from the eventual four-point projection also use the target. None can be supplied to the prediction stage as if independently known.

The frozen three-known-corner benchmark is the appropriate empirical check: save its predictions without target finiteness, component maxima, author means or target-dependent residuals, then evaluate every eligible case. A successful held-out mean-error comparison would support only that declared predictive performance. It would still not prove a global supremum error bound. A failed prediction must not be rescued by replacing it with a post-target projected upper bound, choosing a favorable rank threshold, or calling the observed residual a confidence interval.

## Monotone nuisance transformations

Applying one common strictly increasing map g_k to each **final vertex trait summary** preserves

```text
f_k(v)>=min(f_k(s),f_k(t))
iff g_k(f_k(v))>=min(g_k(f_k(s)),g_k(f_k(t))).
```

Thus the primary endpoint-relative classification is invariant, including its joint intersection. A merely nondecreasing map can merge ties and remove an obstruction; strict monotonicity, or order reflection on these finite values, is necessary. Transforming each preparation before taking arithmetic means is not the same operation and need not preserve the order of the resulting means.

Additivity, projected rank and compromise magnitude are not invariant to nonlinear monotone changes. For example, let both traits initially equal z=x1+2x2. They are additive with rank1. Replace trait2 by the globally strictly increasing map z cubed. In corner order its values are(0,8,1,27), with contrast18. The projected effect columns become(1,10) and(2,17), with determinant-3 and rank2. The primary endpoint-relative classifications are unchanged for every endpoint choice, while exact observed additivity and projected rank change. This alone blocks interpreting interaction or rank as a representation-independent mechanism. Positive affine maps preserve additivity and exact rank, though normalized bound sizes require correspondingly declared scales. Nonlinear maps also change the native-unit nonzero R; the zero-versus-positive accessibility statement remains ordinal.

## Recommended role and stopping decision

Retain the result as a reviewed conditional product-graph consequence and as a synthetic control explaining that opposing additive effects can create joint path conflict. On the present Q2 data, use exact graph classification and the prespecified local prediction benchmark for the empirical claims. If a rank-bound table is later shown, label it an optional post-seal a posteriori description, give its exact scales/domain and residual definition, and compare it to the exact two-route R. It cannot be used as independent confirmation, a latent-noise certificate, a full-space barrier, or a new predictive success.

No current metadata establish a large product of codon states with externally calibrated uniform error. Consequently the headline that this bound supplies a new quantitative prediction of biological compromise on the available sampled components is **not supported**. The mathematical bound remains sound under its reviewed hypotheses; the empirical use is sharply limited. Further generic rank calculations on these tiny products would not resolve the missing biological calibration or graph coverage. Stop this bounded application audit here rather than expanding the theorem or opening outcomes.

`verify_hypotheticals.py` checks the displayed algebra and graph counterexamples exactly. Its finite checks are not a new proof of C5. `manifest.json` records every local source read and the synthetic check receipt. No manuscript, repository or biological outcome file was edited or opened.
