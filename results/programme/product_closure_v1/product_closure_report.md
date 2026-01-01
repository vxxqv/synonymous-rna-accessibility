# Exact additive product queries: current result

The proposed closure survived the bounded falsification study. A direct proof and an exact implementation show how to answer a threshold connectivity query using explicit factor graphs without constructing their Cartesian product. A separate independent proof review and explicit-product checks support the theorem; the witness reconstruction adds quadratic factor-size time beyond the decision query. Generic composition of additive merge trees is established prior work, so a novelty claim is premature.

The correct signature consists of the final coordinate projection sets obtained by allowing each factor to expand while the others sit at reachable peak values. Two retained endpoints are connected exactly when their closures agree. The procedure also constructs a common reachable peak anchor. A factor frontier implementation uses O((N+M) log(max(2,N))) arithmetic and comparison operations per starting endpoint, where N and M are sums of factor sizes. Exact rational widest capacity follows by denominator normalization and integer threshold search, with polynomial bit costs. Explicit paths have length at most twice the sum of (n_i-1)^2 for connected endpoints.

The component itself need not be the retained part of the product of its projections. Two paths scored (2,0,1), at threshold 2, supply the required counterexample: the component of (0,0) projects onto every vertex of both factors, while retained (2,2) is isolated. Its own closure distinguishes it. The implementation rejects the incorrect one-sided projection test.

| Finite domain | Factor and score datasets | Threshold worlds | Retained-start closure queries |
|:---|---:|---:|---:|
| Every binary-scored simple factor on up to three vertices, in unordered pairs | 2,775 | 13,875 | 67,065 |
| Two ternary-scored P3 factors | 378 | 2,646 | 13,608 |
| Three ternary-scored P3 factors, modulo factor permutation and independent reversal | 1,140 | 10,260 | 153,900 |
| Four binary-scored P2 factors, modulo factor permutation | 35 | 245 | 2,240 |
| Total | 4,328 | 27,026 | 236,813 |

For every retained start, the factor algorithm agreed with explicit product BFS on the exact coordinate projections and complete component partition. A separately written cyclic closure also agreed. Every returned peak anchor was canonical, reachable by its recorded path, and maximal over the explicit component. Every path edge and score threshold was checked. The operation counters respected single admission and single adjacency scan bounds. Symmetry reductions use explicit graph automorphisms; finite counts are the representatives actually executed.

An additional 5,184 rational widest-capacity queries, covering every ordered endpoint pair in the stated P3 products, matched an explicit product widest-path oracle. All corresponding endpoint witnesses were valid. Boundary checks covered disconnected factors, a threshold with a new denominator and a denominator exceeding 100 bits. The computation completed in 17.1 seconds, and its manifest validates. These finite assertions support implementation correctness in the recorded domain; the universal statement rests on the separate proof.

The nearest checked prior is [Baryshnikov (2025), Appendix B, Proposition 7.1](https://link.springer.com/article/10.1007/s41468-025-00224-w), which establishes merge-tree composition for additive functions and describes the associated bars. The article delegates detailed attachment rules to earlier work by Le and gives no explicit complexity bound for the individual endpoint query studied here. Search-indexed material from Le's 2023 thesis, Section 5.1, describes construction for finite sums, but its landing page and PDF returned 403 during this audit. Its full algorithm has therefore not been compared.

The possible additional contribution is an explicit query procedure, complexity accounting and constructive witness, rather than the already known compositional principle. That distinction is potentially useful because a full product merge tree can have exponentially many leaves: d paths scored (1,0,1) have 2^d isolated retained vertices at threshold d. Priority for the query procedure remains unresolved until the closest prior construction is examined fully.

The assumptions are consequential. The factors must be undirected and explicit, moves must change one coordinate along a factor edge, and the score must be additive. The paths need only stay above threshold and may decrease in score. Nothing here establishes monotone adaptation, an interaction model, or a connection to measured biological function. No biological data or intervention analysis was performed.

Files: `src/product_closure.py`, `analysis/programme/product_closure_proof.md`, `analysis/programme/product_closure_checks.py`, and `analysis/programme/product_closure_contract_v1.json`. The evidence is `results/checks.json` and `manifest.json` in this result directory. Reproduce in a new directory with:

```text
python analysis/programme/product_closure_checks.py --output work/product_closure_reproduction
python -m pytest tests/test_product_closure.py -q
```
