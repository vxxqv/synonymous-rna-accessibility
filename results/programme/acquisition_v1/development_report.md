# Exact query benchmark: development result

The proposed path/cut alternation does not meet the frozen efficiency criterion in development. Its median query count equals LazySP forward and its mean is higher. The implementation returns sound certificates in every executed case. No validation seed has been evaluated, and no policy has been adjusted in response to these outcomes.

The experiment includes all 760 frozen development worlds: seven graphs, 38 graph and field combinations, ten seeds, and two initial observation conditions. Each world has five policy evaluations. Endpoints are selected without scores and forced present. Initial extra observations are selected without outcomes. Costs count additional internal vertex measurements, including unsuccessful queries.

| Policy | Median queries | Mean queries | Maximum queries | Mean excess over clairvoyant certificate |
|:---|---:|---:|---:|---:|
| Uniform random | 9 | 12.578 | 61 | 9.618 |
| Highest degree | 7 | 12.811 | 62 | 9.851 |
| LazySP forward | 4 | 5.387 | 47 | 2.428 |
| LazySP bisection | 5 | 6.868 | 51 | 3.909 |
| Path/cut alternation | 4 | 5.821 | 41 | 2.862 |

The path/cut policy wins 89, ties 388 and loses 283 comparisons against LazySP forward. Its mean paired excess is 0.434 queries. Its smaller maximum in this finite sample does not establish a smaller worst case over arbitrary worlds. Against LazySP bisection it wins 207, ties 428 and loses 125 comparisons, with a mean saving of 1.047 queries. These comparisons include all worlds rather than a favorable subset.

Failure cases remain in the complete results. On `independent_path_7:bernoulli_0.5:1004:endpoints`, path/cut uses 17 queries and forward uses six. On `Q6:bernoulli_0.25:1003:endpoints`, the respective counts are 22 and 11. Conversely, on `Q6:bernoulli_0.5:1002:quarter`, they are eight and 35. This heterogeneity is a reason to retain the strong forward baseline, not evidence for selecting a favorable graph class after inspection.

All 36,833 simulated query prefixes were checked against complete truth and an independent widest path certificate implementation. There were zero false certificates. Exhaustive reference checks covered all 729 ternary partial observation states on Q3, all 4,096 compatible state and completion pairs, and 351 unresolved states with brute force minimum separator comparisons. Clairvoyant path/cut costs agree with exhaustive certificate enumeration. Exact Q3 minimax recursion gives six additional queries from endpoints alone. That minimax value is a worst case reference; the clairvoyant bound is conditional on knowing the truth and is not an implementable policy.

Five hand-computable unit tests also pass. They check unknown observations, alternative paths, protection of known present vertices, adjacent endpoints, invalid observations, and the distinction between a true certificate and worst case query cost. The bounded experiment finished in 14.7 seconds. The manifest status is `completed`, and the manifest validator accepted all recorded file hashes. These checks validate the stated finite calculation, not a universal policy improvement theorem.

The full graph and field strata are in `development/summary.json`; all truth worlds are in `development/worlds.json`; every policy cost and complete query order is in `development/policy_results.json`. This is abstract synthetic graph evidence. It is neither biological validation nor a model of experimental measurement error. Uncertainty across future validation seeds would describe this fixed simulation design, not a biological population.

The justified decision is to retain the exact certificate software and the negative development result. Validation seeds 2000 through 2049 remain untouched. A validation run requires review of the fixed policies and source hashes, a separate execution record, and no tuning or exclusions based on heldout outcomes. The existing benchmark entry point intentionally executes development only.

Reproduction from the repository root:

```text
python analysis/programme/acquisition_benchmark.py --output work/acquisition_development_reproduction
python -m pytest tests/test_acquisition.py -q
```

The original bounded run and input hashes are recorded in `development/manifest.json`. Query costs and query orders are deterministic under the recorded versions. Runtime fields vary between machines and are not scientific endpoints.
