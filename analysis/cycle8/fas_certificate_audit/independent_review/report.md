# Independent FAS short-certificate review

**Verdict: computationally verified in the stated finite range, with the persistence and E4 algebra proved under the explicit restrictions below. No discrepancy or counterexample to the stated restricted claims was found.** This is a separately authored review with inherited project context, not a blinded external review. The author result was visible before independent implementation. The reviewer did not read, import or execute `analyze.py`.

## Claim and evidence

The objects are the fixed 72 positive empirical envelopes, the eight-column treatment-coded codon design and all strict comparisons `U_i < L_j`. The question is the minimum number of these comparisons whose design differences have a positive dependence. It is not the minimum number of measurements, a selected biological replicate, a global barrier or a calibrated recovery estimate.

`verify.py` verifies the pinned Cycle 6 source hash before decoding. It reconstructs codon states directly from the archived sequences and declared coordinates `(23,38,29,32,41)`, checks their equality with the archived states, builds the treatment coding independently, and generates all comparisons from the rational endpoints. No source counts or preparation-level fields enter the calculation.

| Obligation | Independent result |
|---|---|
| Original source identity | SHA256 `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466` matches. |
| Complete comparison reconstruction | 74 comparisons; all endpoint identities, eight-coordinate differences and edge-loss thresholds match the author output exactly. |
| Every singleton support | 74 checked, zero positive dependencies. |
| Every pair support | 2701 checked, zero positive dependencies. |
| Every triple support | 64824 checked, exactly one positive dependency. No support was skipped on these data. |
| Unique minimum | Comparisons `(16,46),(42,13),(47,46)`, weights `(1/3,1/3,1/3)`. Involved vertices are `13,16,42,46,47`. |
| Cancellation | Exactly zero in all eight coordinates, including the intercept coordinate. |
| Persistence | Exclusive upper bound `t = 9128860950384/8609574801263`, where `t=gamma^2`. |
| Equality | Only the selected edge `(42,13)` becomes a tie at this bound; it is no longer a forced strict order. |
| Global boundary, separate upper argument | Supplied beta `(0,0,-4,3,4,1,1,-1)` independently satisfies all 50 strict orders at this value of t. All 2894 ordered-latent interval-completion inequalities also pass exactly. |
| Known-answer checks | 45 passed: null vectors, rank-one nonminimal triples, a false dependence revealed only by a third coordinate, minimal-support deletion, and rational product factorizations. |

The three difference rows are

```
(0, 1,0,0,0,0, 1, 0)
(0,-1,0,0,0,0,-1, 1)
(0, 0,0,0,0,0, 0,-1).
```

Their sum is zero. If each difference had strictly positive inner product with a latent coefficient vector, summing the three inequalities would give `0>0`. Positive weights therefore contradict the named latent-order model. No omitted-coordinate projection is used.

## Why the enumeration is exhaustive for its claim

A singleton dependency is a zero row. A two-row dependency is an opposite-ray pair and is checked by an exact null vector in one nonzero coordinate followed by verification in every coordinate. For three rows of rank two, any independent pair of coordinate columns gives the one-dimensional left kernel by their cross product; strict common sign and cancellation in all eight columns decide positivity. A rank-three triple cannot have a dependency. A rank-one triple that has a strictly positive dependence necessarily contains either a zero row or an opposite-ray pair, so cannot be minimal. These alternatives cover all triples. Nonminimal supports are filtered only after the smaller supports have been enumerated.

There exists a three-comparison contradiction and there is no one- or two-comparison contradiction among the full 74-row system. Thus three is the **globally minimum comparison support cardinality**. The computation does not enumerate larger minimal positive dependencies, prove uniqueness of all possible certificates, or minimize genotype/measurement cardinality. The reported uniqueness means one minimum-cardinality support among these 74 comparisons; multiplying its weights by a common positive constant is the same dependence.

## Strict persistence and the global boundary

Let `a_i=0` for fixed reference vertex 46 and `a_i=1` otherwise. Its widened interval is `[L_i/gamma^a_i, U_i*gamma^a_i]`. For positive endpoints, comparison `i<j` remains strict exactly when

`U_i^2 * t^(a_i+a_j) < L_j^2`.

At least one endpoint of each present comparison is unanchored. Therefore the loss threshold is `(L_j/U_i)^2` for an anchor comparison and `L_j/U_i` for a comparison between two unanchored coordinates. The selected thresholds are respectively

```
61905156488289/48559156033600
9128860950384/8609574801263
15507410823025/10617848318016.
```

Their minimum is the previously established t*. The certificate proves infeasibility for every `1 <= t < t*`, not at equality. Equality of box endpoints does not force a strict latent order: a common nondecreasing response may have a plateau.

Loss of one certificate alone would not prove global feasibility. Here that separate gap is closed by exact verification of the archived boundary beta. More strongly, for every ordered pair of latent values `q_i <= q_j`, the verifier checks `L_i/gamma^a_i <= U_j*gamma^a_j` using squared rational inequalities. Group the equal latent values, assign each group the largest lower endpoint of all preceding groups, and extend constantly or monotonically between them. Every assigned value lies inside every corresponding upper endpoint by those pairwise checks. This constructs a compatible common nondecreasing response at equality. The nested boxes then remain feasible for larger t. Thus the lower certificate and independently checked upper witness establish the same sharp global feasibility boundary. This statement relies on arbitrary nondecreasing response functions; it is not a continuity, smoothness or strict-monotonicity result.

## E4: what multiplies and what does not

For finite `gamma,eta >= 1`, suppose measurement and recovery allowances are independent coordinatewise symmetric factors, with **the same coordinates anchored at factor one in both allowances**. On a nonanchored coordinate, their multiplier sets satisfy

`[1/gamma,gamma] * [1/eta,eta] = [1/(gamma*eta),gamma*eta]`.

The inclusion from left to right is immediate. Conversely, for any r in the right interval, choose a measurement factor in

`[max(1/gamma,r/eta), min(gamma,r*eta)]`.

This interval is nonempty, and the remaining factor `r/measurement_factor` lies in `[1/eta,eta]`. Apply this construction independently to every nonanchored coordinate and leave anchors fixed. The measurement-expanded box followed by recovery-expanded completion is therefore exactly the anchored box with widening `gamma*eta`. A recovery divisor instead of multiplier gives the same set because the bounds are reciprocal. Consequently every model-feasibility predicate defined on this one-field box depends only on the product. This is set algebra, not a new identification theorem or evidence of a measured mechanism.

The phrase “independent unstructured” is load-bearing. For example, take point observations `(1,4,1)` with first coordinate anchored, measurement factor one and recovery allowance two. Coordinatewise recovery permits corrected values `(1,2,2)`. If the last two coordinates must share one recovery factor, their ratio stays four and they cannot become equal. The product box incorrectly admits that equality under the structured model. The same issue arises for calibrated unequal coordinate bounds, known pairwise restrictions, differing anchor conventions, or factors coupled by a physical model.

Also, the quantifiers here are **one field in one rectangular completion and one persistent factor per genotype**. They do not allow a different persistent factor for every preparation. The two positive fields `(1,2)` and `(1,1/2)` can each be represented by a constant inclusion field with different recovery factors, but not simultaneously by that same constant inclusion field and one shared factor. Empirical boxes can combine endpoints never jointly observed. E4 does not solve this ambiguity or identify actual assay recovery. Retiring the redundant two-axis unstructured plot is justified; replacing it by a structured biological claim is not.

## Scope, dependencies and remaining limits

The dependency chain is: pinned empirical envelopes and categorical states; independent rational comparison reconstruction; complete finite small-support enumeration; the elementary positive-dependence contradiction; exact widening inequalities; an independently verified boundary primal and isotonic completion; and separate elementary E4 box factorization. The theorem of alternatives is established prior mathematics. The protocol cites Boyd and Vandenberghe, *Convex Optimization* (2004), Example 2.21, page 50, but the reviewer did not reopen that external PDF in this bounded task. The contradiction and completeness within supports up to three are derived above rather than delegated to that citation. No priority claim beyond the inspected sources is supported.

This is an exposed exploratory explanation of the existing assay-order incompatibility. Vertex 46 is the normalized fixed WT identity, not a noisy additional biological observation. The smaller five-vertex certificate was selected after outcomes and cannot replace the previously frozen eight-vertex replication target as if prespecified. It does not distinguish inclusion from total RNA, differential recovery, or genotype-by-preparation effects. No calibrated confidence region, independent source validation, physical barrier or causal mechanism has been established.

No failed enumeration or arithmetic attempt occurred in this review. Existing source and author files remained unchanged. Output `receipt.json` preserves every successful short support, exact weights, cancellation, thresholds, boundary checks and hashes. The verifier refuses optimized Python and existing output files.

Reproduction, from a repository/workspace containing these immutable inputs:

```
python -B verify.py --source full_result_v1.json --author-result result.json --robustness robustness.json --output NEW-receipt.json
```

Use the supplied paths or relocated byte-identical inputs. Python standard library only. The strongest safe new sentence is: “A complete audit of all one-, two- and three-comparison supports found one three-comparison minimum certificate; its strict persistence reaches the independently verified anchored feasibility boundary.”
