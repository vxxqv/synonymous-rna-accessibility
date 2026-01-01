# Separately authored review of U1-U3

## Verdict

The mathematical statements U1, U2 and U3 are **proved as written**, under their stated finite graph, connected endpoint and nonempty bounded rational parameter-set hypotheses. The exact implementation agrees with a separately authored rational oracle in all 81 supplied cases. There is no mathematical counterexample in this review.

Two qualifications affect presentation and computational accounting:

1. The fixture named `two_dimensional_corner_failure` does not demonstrate that failure. Its capacity is `|theta1+theta2|` on the square, so its minimum zero is already attained at two mixed corners. Preserve its results but rename its role. A separate genuine two-dimensional challenge, suggested by the author during this review, is checked below.
2. The stated per-candidate arithmetic expression does not include construction of the hyperplane list. In particular the current implementation still performs quadratic pair construction when `d=0`, although there is one candidate and only near-linear graph work in that expression. Include preprocessing in any total implementation bound, or return the single feasible empty vector before constructing score equalities. This does not affect fixed-dimension polynomial complexity or any answer.

The review is separately authored with inherited project context, **not blinded or isolated**. A new isolated reviewer was unavailable because of the agent-thread limit. The reviewer did not author these claims. The proof was reconstructed before consulting `recorded_check/results.json`. The independent oracle does not import the author's solver or verifier. Agreement between authors is not itself proof.

## Exact review point

Repository base: `e9846dca9f74b2c5c34e984eed970e996b9df2bf`, observed clean during the run.

| Input | SHA256 |
| --- | --- |
| `work/uncertainty_cycle3/claims.md` | `db1e36c8860c32e2b790588d238b5863c595b1cf5b03ca7674a4fc1dd914b027` |
| `work/uncertainty_cycle3/affine_polytope.py` | `b068f2b2b1faeb7b33e2aa15c5803966e52acc500b77e0f5d476569a6046a8da` |
| `work/uncertainty_cycle3/protocol.json` | `0edbc70f95f2b15b2056da9a9f185c5eb6805f1f7e6877966291494b402e70ae` |
| `work/uncertainty_cycle3/cases.json` | `a3cc796df4de6459986112bf08c8c61fb1aa8bbc27955d88467d84b77c45afa5` |

All four were unchanged during execution. The frozen independent script and input hashes are in `freeze.json` and `manifest.json` in this review directory. No author source was edited.

## Dependency and obligation matrix

| Claim or obligation | Dependency type | Finding | Remaining scope |
| --- | --- | --- | --- |
| U1: capacity is affine on each pulled-back weak score-order cell | Internal proof from finite max/min expression | Passed. One fixed coordinate realizes each path bottleneck and the outer maximum; tied representatives coincide. | No need to assume distinct scores. |
| U1: barrier is affine on the same cell | Internal proof | Passed. Both endpoints occur among the ordered scores. | Endpoints cannot be omitted from the arrangement. |
| U1: cell extrema are at enumerated ambient intersections | Polyhedral fact, reconstructed in `proof_reconstruction.md` | Passed. A vertex has spanning active normals even when Q is rank deficient. | Requires bounded nonempty cells, as stated. |
| U1: rational and parameter complexity | Determinant bounds and explicit candidate count | Passed for polynomial time at each fixed d. The exponent depends on d. | The current method is not given an FPT bound. Include preprocessing in total runtime claims. |
| U1: zero dimension, singleton, coincident endpoints, negative and tied scores | Definitions plus exact checks | Passed. The API separately reports empty Q and disconnected endpoints. | No finite barrier is assigned when endpoints are incompatible. |
| U2: possible-connectivity hardness | Polynomial many-one reduction from 3SAT | Passed. Undirected backtracking cannot cross a missing layer or reconcile opposite literals at threshold 1. | Clauses may be restricted to nonempty three-literal clauses for the hardness source. |
| U2: possible-connectivity membership in NP | Simple path plus feasible rational LP witness | Passed. Path length is at most n and a polynomial-bit rational witness exists. | Explicit graph and rational input, not a succinct exponential graph. |
| U2: guaranteed-connectivity hardness | Polynomial many-one reduction from UNSAT | Passed. A simple path must select one internally disjoint clause path; Boolean satisfaction supplies a disconnection witness. | Formula evaluations at fractional parameters are not Boolean assignments themselves. |
| U2: universal direction and ties | Internal rounding argument | Passed. Rounding at 1/2 yields a false clause whose entire path has scores at least 1/2. | Uses closed retained thresholds. |
| U2: membership in coNP | Separator plus rational strict-feasibility witness | Passed. Add common positive slack bounded by 1, maximize it over Q, and use a polynomial-bit positive rational optimum. | Separators must be allowed to contain either endpoint. |
| U3: unique minimal vertex set suffices | Internal min/min commutation | Passed for every nonempty compact set of real score vectors. | No independence, convexity or affine assumption needed in this direction. |
| U3: necessity already for one affine parameter | Explicit construction from two incomparable minimal sets | Passed. Guaranteed capacity is zero and fixed capacity is minus one. | A simple-path count cannot replace minimal vertex-set count. |
| U3: chords and direct endpoint edge | Set minimality and direct finite examples | Passed. A direct edge gives the unique set {s,t} even if longer paths exist. | For s=t the unique minimal set is {s}. |
| Generic normal forms and score-order structure | Established external prior | Passed as a prior limitation. | No new generic lattice-polynomial or arrangement algorithm is claimed. |
| Statistical or biological interpretation | Outside the theorem hypotheses | Out of scope. | No calibrated uncertainty set, confidence level or biological mechanism is established. |

The only external complexity leaf is the standard NP completeness of 3SAT, and hence coNP completeness of UNSAT. The constructions and membership certificates are supplied explicitly here. No unproved equality of complexity classes is assumed.

## Decisive arguments

For U1, the important lower-dimensional issue is ambient rank. If the active normals at a purported cell vertex failed to span the ambient space, a nonzero vector orthogonal to all of them would allow a small movement in both directions without violating any inactive inequality. The point would not be a vertex. Equality constraints represented by opposite inequalities therefore supply the missing ambient directions for a line or point. The candidate set can contain extra points without causing an error: every retained candidate is in Q, and it contains at least one optimizer for every cell.

For U2, membership is not inferred from the reductions. A possible-connectivity certificate contains a simple path and rational theta satisfying its vertex inequalities. A failure of guaranteed connectivity contains a separator and strict low-score inequalities. Because Q is bounded, maximizing a slack epsilon subject to `0 <= epsilon <= 1` gives a bounded rational LP. A strictly feasible separator gives a positive optimum, and basic feasible solution determinant bounds give a polynomial-bit witness. This deals with strict rather than closed failure conditions, including arbitrarily small positive slack.

For U3, distinct minimal sets P and Q are incomparable, and their intersection cannot connect. The proposed affine construction gives scores 1 on the intersection and at most `|theta|` elsewhere. P or Q achieves `|theta|`, so the capacity is in fact exactly `|theta|`. Its minimum is zero. At the coordinatewise lower field every vertex outside the disconnected intersection has score minus one, giving fixed capacity minus one. The construction preserves endpoint scores at one. This proves the quantifier gap without sampling parameter values.

The phrase that extrema occur in T should mean that extremal **values are attained** in T. It does not mean that every optimizer belongs to T. Constant-score examples have a continuum of optimizers. This is a wording clarification, not a failure of the algorithmic conclusion.

## Independent exact computation

The independent implementation enumerates minimal connecting vertex sets and minimal vertex separators directly from vertex subsets. It then uses exact Fourier-Motzkin elimination to solve path hypograph and separator epigraph LPs. It never enumerates score-order cells. The independent barrier LP identities and their derivation were recorded before the author results were opened in `proof_reconstruction.md`.

The initial independent run used Python 3.12.14 and exact `fractions.Fraction`, completed in 0.608 seconds, and made 1,127 LP projections and 2,219 elimination steps. The largest intermediate system contained 19 inequalities. No tolerance, floating reconstruction or randomized oracle was used. The 180-second and 20,000-inequality limits were not reached.

Results:

* All 81 supplied case statuses agree with the author receipt. The 79 connected cases give 395 exact rational comparisons: fixed, guaranteed, possible, minimum barrier and maximum barrier. All agree.
* All 1,098 labelled simple graphs on two through five vertices were checked. There are 928 with connected distinguished endpoints: 810 have one minimal vertex set and 118 have more than one. The latter all exhibit the proposed exact quantifier gap. The oracle also checked 28,780 retained vertex subsets against the minimal-set representation. A singleton graph verifies the coincident-endpoint convention.
* The 255 nonempty subsets of the eight full three-variable clauses were checked, together with four explicitly listed unit, repeated-literal and tautological examples. There are 259 formulas, four unsatisfiable, 2,050 Boolean scenario checks and 6,903 ternary checks including exact 1/2 ties. These corroborate the reductions; the proof handles the continuum and arbitrary formula size.

The initial mathematical receipt was saved before the author receipt was read. `comparison.json` records the subsequent exact comparison and both receipt hashes. The author's original comparison used a floating LP tolerance; this review removes that tolerance for the stated finite cases, without claiming that finite agreement proves U1.

### Additional corner challenge

The author suggested a separate six-vertex graph during review, after the independent initial receipt. Its four parallel intermediate vertices have scores `theta1`, `-theta1`, `theta2`, `-theta2`, endpoints have score 1, and Q is `[-1,1]^2`. This is explicitly an author-suggested additional case, not a blinded reviewer discovery and not part of the frozen 81 cases.

Here `C=max(|theta1|,|theta2|)` and `B=1-C`. Every square corner has capacity 1 and barrier 0; the origin has capacity 0 and barrier 1. The exact LP oracle verifies fixed capacity minus one, guaranteed zero, possible one, and barrier range [0,1]. Both latent coordinates matter. It is a valid replacement for the original fixture's advertised role, while preserving the original result.

## Prior and contribution boundary

[Marichal, *Weighted lattice polynomials*, arXiv:0706.0570v2](https://arxiv.org/pdf/0706.0570), Proposition 2, gives the max-of-mins and dual min-of-maxes forms. Theorem 15 expresses these functions on weak coordinate-order regions. This is exactly the generic structure behind U1 and the minimal-set representation used in U3. The paper works on bounded distributive lattices; compact real score ranges can be placed in one bounded real interval. Its statements do not themselves formulate these graph-specific uncertainty quantifiers or the SAT constructions. Nevertheless, the normal-form machinery and score-order argument are established prior, and U1 is a direct application. U3 is a short structural consequence, not evidence for a new general lattice theory. U2's encodings are explicit and correct, but exact priority against correlated-network optimization literature remains unresolved.

The strongest defensible contribution is a precise accessibility formulation: exact ranges for low latent dimension, a variable-dimension computational boundary under an explicit shared-parameter representation, and an exact criterion for when route adjustment can matter. These results do not identify a statistical uncertainty region or demonstrate a new biological effect.

## Reproduction and closure

From the workspace root, in a fresh copy of this review output directory so existing receipts are not overwritten:

```text
work/env/Scripts/python.exe work/uncertainty_cycle3/independent_review/check_exact.py
work/env/Scripts/python.exe work/uncertainty_cycle3/independent_review/compare_receipts.py
```

`manifest.json` passed the installed Mathbox v2 validator with input and output hash checks. The scripts refuse to replace the existing principal receipts. This review is complete. Before integrating, rename the misleading fixture, add the separately labelled genuine corner example if desired, and clarify preprocessing or implement the d=0 shortcut in a new revision. Retain the old freeze and receipts. No further computation is needed to resolve the reviewed mathematical claims.
