# C5 joint-function mathematical handoff

This is one bounded proof/counterexample package, not a manuscript addition or empirical analysis. The repository baseline is `a2cd60afb96be6bd3b11ffc1c692471f979f3251` and was read only. All new files are in this directory. No biological outcomes were read.

Read `prior_map.md` first. Its principal conclusions are that exact constant-trait joint reachability is already hard, that the generic vector-ordering argument is classical, and that Kirsch2026 is a close conceptual prior rather than an unrelated overlapping-gene paper.

`claim_and_proof.md` states the exact graph, scale, interval and endpoint conventions, then proves C5.1-C5.3. The candidate is a formal product/coding consequence of Steinitz: a rank-dependent compromise bound plus an essential local codon-path excursion. The independent review is pending. Do not call this a new generic vector-balancing or bottleneck theorem.

`oracle_examples.json` supplies ten fully specified finite graphs, exact rational trait vectors, thresholds and expected R values for a separately authored oracle. It includes baseline disconnection, local Leu gates, antagonistic additive traits, background detours, coupled-pairing failure, strict/weak plateaus, interval constraints, coincident endpoints and the sharp2epsilon error example.

`verify_finite.py` is an author implementation with two different path algorithms: exact critical-threshold BFS and minmax Floyd. It also exhausts all729ordered lists of three2D vectors in {-1,0,1}^2, trying all six orders and computing exact distances to the endpoint segment. These checks support implementation and examples in this stated range; they do not prove the universal result. Run with Python3.12.14 and standard-library arithmetic only, into a fresh directory:

```text
python verify_finite.py --out NEW_DIRECTORY
```

The retained run is `finite_run_v1/`, with source hashes before/after and its exact output hash in `manifest.json`. No seed, numerical tolerance, hidden sampling, empirical data, Git dependency or network call is used by this finite verification. Disabled assertions and an existing output directory are refused.

The exact original Steinitz source was not fully retrieved. The load-bearing zero-sum statement was verified in Barany2024 Theorem1.2, a primary research manuscript, which attributes the result to Grinberg-Sevastyanov1980. Other inaccessible sources are explicitly identified in the prior map. No fetched literature PDF or copied full text is included.

Next action: root should send the precise statements/proof and raw examples to a fresh adversarial reviewer before integration. In particular challenge the intrinsic-subspace norm, collinear dummy removal with mixed total signs, local endpoint-box lift, rank-zero case, interval margins, product assumption and the global scope required by2epsilon. No wider theorem development or empirical execution is needed to review this package.
