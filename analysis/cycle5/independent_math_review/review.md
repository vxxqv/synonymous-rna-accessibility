# Independent C5 proof review

**Primary verdict: proved as written, under the stated product/additive/closed-box hypotheses and the explicit inherited hypotheses of the approximation corollary.** The internal argument is valid using the authenticated classical Steinitz theorem. No falsifying counterexample was found under those hypotheses. Three standalone-statement clarifications below should be made before manuscript integration; they do not change the bounds.

This is a separately authored adversarial review with inherited project context, not a blinded or external referee review. Agreement and the finite checks are not the proof. The source argument and the reconstruction below supply the mathematical evidence.

Reviewed revisions:

- `claim_and_proof.md`: `d5ef26923060fd587bb35b0c2d03d6ae39c34fda1e9a967f87e723e00d5cadb1`.
- `prior_map.md`: `252e2ba00b776e5feabe24f91a074dbe09100bf99430a3d548b593d8df9e1b30`.
- Exact fixture and source hashes are in `checks.json` and `manifest.json`.

## Claim cards

| Claim | Quantified objects and conclusion | Verdict |
|---|---|---|
| C5.1, lines27-46 | Finite connected Cartesian factors, real additive scaled traits, nonempty closed axis-aligned acceptable box, eligible endpoints. Endpoint-effect rank r, largest endpoint effect M and minimum local box excursion eta give R<=rM+eta. | Proved as written. |
| C5.2, lines48-50 | Same objects; each endpoint has at least sigma inward slack for each finite coordinate bound. Then R<=max(0,rM+eta-sigma). | Proved as written with the stated inward-coordinate meaning of margin. Spell it out algebraically. |
| C5.3, lines82-86 | Same product graph and box; additive model is uniformly within epsilon of the true scaled traits at every graph vertex; true endpoints eligible. Then R_true<=rM+eta+2epsilon, with r,M,eta computed from the additive model. | Proved as a corollary of C5.1. Explicitly restate inherited product hypotheses if quoted separately. |

The bounds concern the declared graph and scaled trait requirements. They do not identify biological validity of additivity, error calibration, trait scales, assay thresholds, complete coverage, evolutionary fitness or mechanism.

## Dependency graph

1. The violation function equals infinity-norm distance from the scaled trait vector to the closed acceptable box.
2. Product domain and additivity identify each completed factor change with one endpoint-effect vector.
3. Classical Steinitz in the intrinsic endpoint-effect span gives a bounded zero-sum ordering after dummy augmentation.
4. Removing dummies gives boundary proximity to one common endpoint segment.
5. A rectangular acceptable set plus a local endpoint-box excursion bound controls the full nucleotide lift.
6. Inward coordinate margin sharpens the distance bound.
7. Uniform model error expands the endpoint box once and the resulting true route once more.

Steps1,2,4-7 are internal arguments checked below. Step3 is the external theorem, verified in Barany's primary author manuscript. The computational fixtures are downstream corroboration, not a leaf used to prove a universal implication. No circular dependency was found.

## Obligation matrix

| Obligation | Result | Evidence |
|---|---|---|
| Intrinsic rank rather than ambient trait count | Passed | Restriction of the ambient infinity norm to the real span V is an r-dimensional norm with a centrally symmetric relative unit ball. |
| Dummy norm and membership | Passed | d=sum(v_i) belongs to V; ||-d/n||<=M. Labels preserve repeated vectors. |
| Common interpolation parameter after deletion | Passed | Each original prefix has one integer dummy count j, hence one scalar j/n for all coordinates. |
| Local paths may leave V | Passed | Only block boundaries use the rank argument. Local interiors are controlled in ambient coordinates by eta, not assumed to remain in V. |
| Endpoint box lifting | Passed | Coordinatewise bounds from the two block boundaries suffice because K is a product of intervals. |
| n=0, M=0 and rank0 | Passed | n=0 gives a constant route. M=0 gives identical boundary phenotypes; eta remains mandatory. No division by n or norm rescaling is used in the exceptional branch. |
| Margin and unbounded coordinate sides | Passed | Apply sigma only to existing finite bounds. Equality-width coordinates force sigma0. |
| Approximation and sharp2epsilon | Passed | Both endpoint and interior errors are required; the Leu construction attains2epsilon for epsilon>0. |
| Nonproduct compatibility constraints | Passed as exclusion | The independent-set counterexample violates the product assumption and refutes extending the bound merely from additive scores. |
| Global optimum versus block witness | Passed | The supplied Q3 detour has zero compromise even though every endpoint-block order incurs one. |
| Strict thresholds | Passed as separate convention | A positive-part violation cannot encode strict connectivity at a plateau; signed violation and an open critical-value condition are needed. |
| Exact versus finite verification | Passed | All ten fixtures were recomputed independently; the report does not treat them as the general proof. |

## Reconstruction of the decisive implications

### C5.1

For n>0 and M>0, put V=span(v_i) and d=sum(v_i). Each of the n dummy vectors -d/n lies in V and in the norm ball of radius M. The2n-vector list sums to zero. Steinitz applied after division by M yields partial sums of norm at most rM in V. No conversion through a Euclidean norm or change-of-basis distortion is needed.

Consider the augmented prefix ending exactly at the kth original vector. If j dummies preceded it, then its sum is `S_k-(j/n)d`, where S_k is the original-vector prefix. Consequently

`F(s)+S_k = F(s)+(j/n)(F(t)-F(s))+e_k`, with `||e_k||_infinity<=rM`.

The scalar j/n is common to the entire vector and lies in[0,1]. It is not necessary that j/n=k/n, that the endpoint displacement be componentwise positive, or that all dummies lie at the end. Both the empty boundary and the final endpoint are already eligible. Every block boundary therefore lies in the rM expansion of K.

During block i, all other factors are frozen. Let its before/after global phenotypes be A and B. A local state z has, for each coordinate k, a value in

`[min(A_k,B_k)-eta_i, max(A_k,B_k)+eta_i]`.

This follows directly from distance to the local endpoint box, with the contribution of the other factors added to all three values. Because A and B lie in K expanded by rM, their entire coordinatewise box lies there too. The intermediate is thus in K expanded by rM+eta. This proves the claimed upper bound. The argument would not follow from convexity alone for a general nonrectangular K.

The explicit local-box example `(0,0) -> (0,1) -> (1,1)` has eta_box0 while its middle point is distance1/2 from the endpoint segment and lies outside the endpoint-effect line. The proof correctly avoids assuming segment proximity for factor interiors.

### C5.2

Define the inward box using `L_k/a_k+sigma` on finite lower bounds and `U_k/a_k-sigma` on finite upper bounds. The endpoint segment stays in this box. The boundary and local-path arguments give inward slack at least `sigma-(rM+eta)`. Violation is the positive part of its negative, exactly `max(0,rM+eta-sigma)`. Unbounded sides impose no extra condition. If sigma exceeds the route error radius, the conclusion is zero violation, not a negative R.

### C5.3

For every vertex x, `||F_true(x)-F_add(x)||_infinity<=epsilon`. Eligible true endpoints imply model endpoints in K expanded by epsilon. C5.1 applied to this expanded box gives a model route within a further rM+eta. Returning to the true values costs another epsilon. Coordinate expansions add, giving the stated rM+eta+2epsilon.

The two epsilon terms cannot generally be merged. With the model constantly -epsilon and the true Leu endpoints0, all other Leu codons -2epsilon are still uniformly epsilon from the model. UUA and CUC are not adjacent, so every route traverses a state with violation2epsilon, even though model r=M=eta=0. Controlling errors only at endpoints would not protect any interior state. A training RMSE, average held-out error or singular-value residual is not the uniform hypothesis.

## Boundary and computation evidence

`check_examples.py` uses exact rational min/max Floyd updates, initialized with both endpoint vertex costs on each edge. It then independently checks the critical value with threshold graph traversal. The author implementation was neither imported nor executed.

| Supplied example | Independently computed R |
|---|---:|
| Leu local excursion at rank0 |3|
| Ser different base components |infinity|
| Rank1 conflict |1|
| Rational effects |1/3|
| Background detour |0|
| Coupled independent-set domain |2|
| Weak plateau |0|
| Equality interval |1|
| Coincident endpoints |0|
| Sharp uniform-error example |1/2|

Every listed edge set also matches the complete induced Hamming-one graph on its supplied labels. This does not itself verify a biological clone or translated peptide. Four additional rational margin cases give R=1,3/4,0,0 for sigma=0,1/4,1,2, respectively. The local-box example and the positive2epsilon example were checked separately. No discrepancy occurred.

Command: `work/env/Scripts/python.exe -I -B work/joint_function_cycle5/independent_math_review/check_examples.py`. Both runs passed. The second run added an exact identity assertion connecting the independently constructed2epsilon values to the supplied fixture; no source theorem or numeric tolerance changed.

## Source, complexity and priority checks

- [Barany2024 author manuscript](https://real.mtak.hu/212570/1/2308.10102v2.pdf), Theorem1.2, PDFpage2, supplies the required zero-sum arbitrary-norm bound d. Its introduction attributes that constant to Grinberg-Sevastyanov. The norm/rank translation used here matches its hypotheses. I did not independently reread the1980 original-language proof.
- [Cardinal et al., section4, Theorem4.7](https://arxiv.org/html/1805.04055v1) explicitly proves constant-inequality constrained-cube reachability PSPACE-complete. Its displayed construction uses23*3+1=70 tuple dimensions, with46 lower and24 upper requirements. The AAA/AAG factor is exactly a binary edge, so repeated Lys factors reproduce the cube without extra genotype vertices. This is a direct notation transfer, not new generic hardness. Fixed rank at most70 does not contradict an existential relaxation bound.
- [Ito-Demaine erratum](https://erikdemaine.org/papers/SubsetSumReconfiguration_JOCO/erratum.pdf) replaces the strong-NP-hardness claim with NP-hardness and gives the corresponding maxmin corollary. The prior map correctly retains this correction.
- [Ambrus-Heck2026](https://doi.org/10.1112/mtk.70085) explicitly studies non-zero-sum ordering and its relation to Steinitz constants. The publisher source reader exposed its framework; a separate local GET403 is preserved, not bypassed. This review does not authenticate every displayed formula or establish an exact priority comparison with its non-zero-sum refinements.

Explicit-graph bottleneck optimization is polynomial in the explicit graph and rational input encoding. That graph may be exponential in coding length. The proof does not supply or claim a new efficient Steinitz ordering algorithm; n! enumeration is not polynomial. Exact rational rank/M and six-vertex codon eta queries are consistent with the stated local computations. No hidden full-genotype enumeration is used in the existential proof itself. Biological application still needs valid contribution tables and global error certification.

This review did not freshly re-audit every source mentioned in `prior_map.md`, including Kirsch, Mohri or the older graph-convexity literature. The report's generic bottleneck, induced-Laplacian and nonmedian cautions are mathematically consistent. No global priority, first use, optimal higher-rank constant or new general algorithm is established. The strongest safe contribution label remains a coding/product-graph corollary of established vector balancing, with an explicit local-intermediate correction and conservative deterministic approximation allowance.

## Precise integration clarifications

1. **Required when C5.3 is quoted alone, line82:** begin with “Under the same Cartesian-product, connected-factor and rectangular-acceptance hypotheses as C5.1”. Additivity on an arbitrary constrained graph is insufficient. With epsilon0, the K2,2 independent-set example would otherwise give the false2<=1 conclusion.
2. **Required notation clarification for C5.2, line48:** define margin by the coordinate inequalities above, omitting infinite bounds. “Distance from every finite face” must not be interpreted as distance only from bounded faces or box vertices.
3. **Required standalone calibration clarification, lines82-86:** replace “every relevant base-graph vertex” with “every vertex of the declared product graph” or explicitly name the complete endpoint component. Do not permit a measured subset to define relevance after seeing outcomes.
4. **Scientific wording at line80:** use “a certified full-graph lower bound on R exceeding this valid upper bound” in place of a potentially sampled-graph “empirical obstruction”. Even then it rejects the conjunction of model, graph and measurement assumptions; it does not identify the restrictive mechanism.

Items1-3 make inherited/implicit meanings explicit; they are not new restrictions on the proof reviewed here. No theorem-bound correction or fixture correction is required. If any of these hypotheses is weakened, this approval does not transfer.

## Remaining limits and stop decision

The decisive general argument passed. The finite tests are limited to the ten archived examples and six declared boundary diagnostics. The external source check verifies the exact Steinitz leaf and nearby hardness/correction boundaries, not exhaustive novelty. There is no justification yet for treating an empirical additive fit, low rank approximation, sampled obstruction or nominal read precision as a calibrated biological instance of these hypotheses. Stop this bounded review here; further work should target a genuine data/model identification gate rather than additional generic proof enumeration.
