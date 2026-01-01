# FAS ordinal and interval-geometry implementation review

**Production numerical verdict: PASS on corrected `ordinal_geometry.py` SHA256 `8872d6452afe4ccbf84e61450d4b1bfbdadf26469f20c0e24452f3abc205f9ac`, within the stated finite integer-design and valid-interval contract.** A genuine false barrier caused by consuming an edge iterator was found and fixed before this verdict. The original failing source and reproducer are retained.

**Test-evidence action:** the reviewed `verify_synthetic.py` SHA256 `85a1b2fdba3a56d526cf941db9d7aa2e6d1c91b7b80732b8e9998c0228bdb4fb` uses removable `assert` statements and does not refuse `-O`/`-OO`. Its normal-mode16-check run passed in an isolated copy. This does not establish that an optimized invocation meaningfully runs those assertions. Add an explicit pre-import optimization refusal before treating arbitrary invocations as reproducibility evidence, or strictly constrain and verify normal-mode execution. No biological outcome was used in this review.

## Independent evidence

The production files were copied into `candidate_v1/` and hashed. No shared source or result was edited. The copy happened after the iterator correction, so its exact hash, not the directory suffix, identifies the reviewed version. The original `e94fba1d...` version was reconstructed by reversing only the two added edge-materialization lines; its complete SHA256 matches the pre-correction source observed in the minimal reproducer. It is retained in `history_iterator_bug/`.

The independent oracle does not import the theorem agent's numeric implementation. In two dimensions it detects the origin in the convex hull of forced-order difference rows by exact rational supports of size one, two or three. Caratheodory's theorem and the separating-hyperplane alternative justify this finite oracle. All returned feasible responses are independently checked against the original intervals and latent features. All returned infeasibility certificates are independently checked for valid forced-pair indices, the original difference vectors, nonnegative rational weights, unit weight sum and zero weighted difference.

| Check | Scope and result |
|---|---|
| Independent ordinal oracle |1,107 cases: four rank/duplicate-feature configurations with216 interval assignments each, plus81 Q2 point fields under three row permutations. All agree:937 feasible and170 infeasible. Unknown upper bounds and touching closed intervals are included.|
| Full72 positive fields |Three row orders of an explicitly additive latent field passed through a nonlinear nondecreasing plateau. Exact returned responses and intervals checked.|
| Full72 negative fields |Three row orders of a four-vertex reciprocal-order obstruction with68 unknown vertices. All return independently valid rational duals.|
| Additional degeneracy/scaling |Repeated and linearly dependent columns, duplicate rows, constant responses, and rational response scales from10^-80 to10^80. No status or certificate mismatch.|
| Widest paths |31,104 unordered endpoint checks against threshold BFS: all64 labelled graphs on four vertices, including disconnected graphs, and81 fields per graph.|
| Interval bounds |6,144 endpoint checks against independent lower/upper threshold connectivity:256 interval assignments on four graphs, including unbounded uppers. Exact positive, zero, disconnected and tied boundaries agree.|
| Edge representations |Lists, reversed lists, duplicate edges, self-loops and one-use iterators give the checked unknown-bridge result after correction.|
| Invalid inputs |Noninteger features, reversed intervals, invalid edge indices and NaNs are refused.|
| Supplied suite |All16 supplied checks pass under normal Python in the isolated review copy.|

`independent_v1.json` holds the main counters; `targeted_results.json` holds19 additional cases and the full72 negative certificates. The tests use exact rational comparisons wherever a certificate is asserted. No random subset, biological data, model fit or outcome-dependent selection is involved.

## Resolved false certificate

On the path0-1-2, take lower bounds `[2,0,2]` and upper bounds `[2,3,2]`. The optimistic endpoint bottleneck is2, so endpoints0/2 have no uniformly excluded accepted-threshold interval. Before correction:

```
bounded_pairs([2,0,2], [2,3,2], iter([(0,1),(1,2)]))
```

returned upper bottleneck `-inf`, `certified_nonempty=True` and width `inf`. The validation loop exhausted the edge iterator; subsequent graph passes saw an empty graph. The identical list input returned the correct upper2, false nonemptiness and width0. Materializing `edges=list(edges)` in both `widest_paths` and `bounded_pairs` resolves both the validation reuse and the lower/upper reuse. The corrected reproducer and the broader independent checks pass. This bug did not affect the original `graph_edges` list workflow, but it was a real unsafe API behavior, not numerical tolerance.

## Why the returned certificates are sound

For closed finite-lower intervals and optional infinite upper bounds, an obligatory order occurs when `U_i<L_j`. A nondecreasing response then requires latent `z_i<z_j`. With finitely many constraints, a strictly feasible latent vector can be scaled so every difference margin is at least1. The code checks its reconstructed rational vector against every original difference before asserting feasibility. The prefix-lower construction gives a single response value for each latent tie and checks membership in each interval. It allows a nondecreasing response with plateaus, not only a strictly increasing response.

An infeasibility return includes nonnegative unit-sum rational weights whose weighted forced-order difference is exactly zero. Multiplying positive latent margins by those weights would give both a positive value and zero, a contradiction. Thus an asserted return is sound even though a floating optimizer proposes the support. A missing or numerically unrecoverable support causes an error rather than a scientific infeasibility label. Rational reconstruction with a finite denominator cap can fail on a feasible problem; this implementation is not claimed to be a complete exact-arithmetic solver for all arbitrarily conditioned integer designs.

For any admissible score vector `L<=f<=U`, max-min path capacity is monotone, giving `C(L)<=C(f)<=C(U)`. If `ell=min(L_s,L_t)`, both endpoints pass at every threshold at most `ell`, while they must be disconnected whenever the threshold exceeds `C(U)`. The correct certified set is therefore **`(C(U),ell]`**, with strict nonemptiness `C(U)<ell`. Equality is not a barrier. The reported width is the width of this guaranteed threshold interval, not a confidence-interval width or the full sharp range of possible endpoint barriers. An unknown upper bridge is allowed arbitrarily high and correctly prevents unjustified exclusion.

## Scope and residual contracts

- The categorical graph helpers assume state tuples of common arity and valid declared factor labels/edges. The current fixed72-state construction meets this condition. `features`/`graph_edges` do not independently authenticate a biological sequence map; the strict reader and frozen metadata must do that. In particular, Python `zip` truncates malformed unequal-length states, so such malformed states must not reach these helpers.
- Integer feature inputs and finite rational lower bounds are required. `None` is the designated infinite upper bound. A production optimization or rational reconstruction error must remain an unresolved computational gate, never be relabelled a biological rejection.
- Disconnected underlying graphs can yield infinite threshold widths; this is a topological disconnection, not an infinite biological effect. The declared72-state graph is connected.
- Ties, finite interval overlap and component uncertainty concern stipulated deterministic score sets. Nothing in these calculations calibrates experimental uncertainty or identifies a mechanism.
- The source helper and source protocol were written by another lane. This is independent implementation review with inherited project context, not a blinded or external review.

No additional numerical defect was found in the corrected source within the stated contract. The original failure is preserved, rather than replaced by an unqualified first-pass success.
