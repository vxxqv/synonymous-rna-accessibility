# Heldout query benchmark

The frozen efficiency criterion failed. Path/cut alternation required a median of five additional queries, compared with four for the forward selector. Its mean paired excess over forward was 0.398 queries, with an approximate 95% interval of 0.261 to 0.535 across the 50 complete simulation seed blocks. The exact certificate implementation passed all executed safety checks, but the proposed policy does not establish an improvement over the strongest included baseline.

All 3,800 heldout worlds and 22,800 policy runs were retained. They cover seven graphs, 38 graph and field combinations, two initial observation conditions and seeds 2000 through 2049. There are 2,652 connected and 1,148 separated worlds. Each policy sees the same world and initial observations. No policy, endpoint, exclusion or success criterion was changed after heldout evaluation began.

| Policy | Median queries | Mean queries | Maximum observed queries | Mean excess over clairvoyant certificate |
|:---|---:|---:|---:|---:|
| Uniform random | 9 | 12.831 | 62 | 9.858 |
| Highest degree | 8 | 13.079 | 62 | 10.106 |
| Forward, vertex adaptation | 4 | 5.420 | 51 | 2.447 |
| Path midpoint | 5 | 6.686 | 49 | 3.713 |
| Path/cut alternation | 5 | 5.818 | 42 | 2.845 |
| Canonical bisection, vertex adaptation | 5 | 6.277 | 48 | 3.304 |

The original five policies and their success criterion were retained. Primary source review corrected the historical `lazy_bisection` reporting label to **Path midpoint**. Canonical bisection was added as a separate sixth baseline before validation, and evaluated first on the previously recorded development worlds. Its vertex adaptation selects the unknown path position farthest from the nearest evaluated position, following [Dellin and Srinivasa, Algorithm 2](https://www.ri.cmu.edu/pub_files/2016/3/paper-lazysp.pdf). The original receipt is preserved. This amendment was made after initial development results and is documented as such.

| Comparison, path/cut minus baseline | Mean difference | Approximate 95% interval | Wins | Ties | Losses |
|:---|---:|:---|---:|---:|---:|
| Uniform random | -7.013 | -7.213 to -6.814 | 3,059 | 426 | 315 |
| Highest degree | -7.261 | -7.463 to -7.059 | 3,059 | 671 | 70 |
| Forward | 0.398 | 0.261 to 0.535 | 470 | 1,963 | 1,367 |
| Path midpoint | -0.868 | -0.995 to -0.742 | 974 | 2,141 | 685 |
| Canonical bisection | -0.459 | -0.571 to -0.348 | 900 | 2,165 | 735 |

Negative differences favor path/cut. Intervals concern the mean paired difference only. Each seed block contains the same complete mix of 76 graph, field and initial observation worlds. Student t intervals use the 50 block means and 49 degrees of freedom. They estimate Monte Carlo uncertainty for this fixed simulation design. They do not quantify biological uncertainty, variation over an unspecified population of graphs, simultaneous coverage across policies, or uncertainty in medians and extremes. An independent calculation using the recorded costs and a separately evaluated Student t quantile reproduced the forward comparison interval.

Failures are retained and can be consequential. On `Q6:bernoulli_0.25:2046:quarter`, a separated world, path/cut used 30 queries and forward used six. On the connected world `Q6:bernoulli_0.5:2010:endpoints`, the respective counts were 37 and 17. The lower maximum observed cost for path/cut does not establish a smaller worst case on arbitrary inputs. Selected successes do not change the failed primary criterion.

All 213,220 query prefixes were checked against complete truth and a separate widest path implementation. No false certificate occurred. Exhaustive checks again covered 729 Q3 partial states, 4,096 compatible state and completion pairs, and 351 brute force separator comparisons. The exact Q3 minimax reference is six queries from endpoints alone. Clairvoyant certificate sizes are lower bounds conditional on knowing the complete truth, not implementable predictions.

The experiment uses synthetic occupancy fields. Forcing both endpoints present is an occupancy edit. In particular, an additive-derived field after this edit is not asserted to remain an exact superlevel restriction of an additive score. This distinction prevents artificial endpoint edits from being interpreted as counterexamples to additive connectivity theory. No biological measurements were used.

The input files matched committed revision `45bc296d84ade7efa824fa0b48a5e9dcefa57b1f` before validation. The freeze SHA256 is `8ec8ce9602bfeb53c7c692854fade45dce8f053168a5d860025fa2e9e1122898`. The bounded run completed in 80.5 seconds, and the manifest validator accepted all declared input and output hashes. The manifest also records the full repository as dirty; all ten declared analysis inputs were separately verified against that exact committed revision before execution. Nine targeted unit tests had passed before the freeze, including output overwrite rejection, hash mismatch rejection, the difference between canonical bisection and midpoint, and seed block uncertainty.

The defensible contribution is an audited certificate workflow and a heldout negative comparison. The result argues for retaining the established forward baseline. It does not support presenting path/cut alternation as a general efficiency advance. No further policy tuning or expensive biological application is justified by this comparison alone.

The complete results are in `validation/results/summary.json`, `worlds.json`, `policy_results.json` and `graphs.json`. Every cost, complete query order, simulated truth world, stratum and seed block mean is retained. The manifest is `validation/manifest.json`.

To reproduce the scientific results from the frozen checkout, use a new output directory:

```text
python analysis/programme/acquisition_benchmark.py --phase validation --output work/acquisition_validation_reproduction
python -m pytest tests/test_acquisition.py -q
```

The driver rejects changed frozen inputs or an existing output directory. Runtime values may differ by machine. Query counts and query orders are the scientific outputs.
