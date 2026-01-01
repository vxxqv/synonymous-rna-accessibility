# D1 exploratory nested-model result

The fixed exploratory test completed once after its implementation freeze. All 60 joint cases and 80 individual-context cases have exact rational certificates. The nested additive model does not exclude any additional case beyond the previously exposed unrestricted common-order result. Context-specific additivity remains compatible in every declared comparison.

| Widening factor | Shared-additive exclusions | Context-specific additive exclusions | Saturated pair-reference exclusions | Prior unrestricted common-order exclusions |
|---|---|---|---|---|
| 1 | Exons 17 and 27 | None | None | Exons 17 and 27 |
| 5/4 | None | None | None | None |
| 3/2 | None | None | None | None |
| 2 | None | None | None | None |

The complete output retains 138 feasible witnesses and two normalized rational contradiction certificates. All 132 state/cell intervals were retained. Every model and widening comparison, the expected saturation checks, the independent-context conjunction and the implication to the weaker common-order reference passed. Increasing the widening factor never changed a feasible model into an infeasible one.

The two original-interval exclusions express incompatible common latent order across the source cell contexts. They were already identified by the weaker B1 comparison and are not a new consequence of additivity. Their disappearance at the first widening factor must accompany any report of them. A widening factor is a hypothetical descriptive sensitivity, not an estimated measurement correction or confidence level.

The sole structurally informative within-context additive constraint lies in exon 17. It is compatible with the intervals in both cells at all four widening levels. This is non-exclusion, not support for additivity. Exons 25, 27, 39 and 42 have full-row-rank descriptor designs within each context, so their context-specific model is incapable of excluding an arbitrary finite order. The saturated observed-pair model is likewise a structural reference; its feasible fits do not identify individual interactions.

The conclusion is bounded: this selected local-state design and these descriptive intervals do not furnish a new interaction exclusion. No complete haplotypes, protein mechanism, independent assay replication or causal marker effects are established. The earlier finite marker-context results and all prior failures remain unchanged. No source values, groups, model choices or widening factors were changed after this run.

## Evidence and review status

The immutable pre-fit freeze is `freeze.json`, SHA-256 `7b0492085b009303790fb1994d02fcfe8bf0760f48dd7eace12a78b3fac1c13e`. It pins the reviewed protocol, fixed matrices, exact B1 result, implementation, synthetic receipt and runtime. The unchanged source runner is SHA-256 `2a58c52ee80ea234bc1466a207ae5e4f9b020f142088721b239abcf1a5b86cc9`.

`actual_v1/result.json` contains all cases, coefficients, latent scores, response completions, forced-pair margins and duals. `actual_v1/receipt.json` pins the result and records the execution. `actual_v1/invocation.json` records all supplied paths and the authorization hash. `actual_v1/progress.json` preserves the completed-case checkpoint. Source outcome fields outside the declared interval contract were skipped without value decoding.

The 566 pre-fit synthetic assertions include independent finite gap-cone comparisons, but the production exact certificate checks are part of this implementation lane. Root's separately authored reconstruction is still required before admission into the manuscript or computational archive. This report does not claim that review has already occurred.

## Replay command

From the workspace root, choose a fresh output directory:

```
python work/research_cycle8/context_design_audit/nested_model/runner.py --freeze work/research_cycle8/context_design_audit/nested_model/freeze.json --authorization work/research_cycle8/context_design_audit/nested_model/authorization.json --b1 work/research_cycle8/cprime_results/primary_v1/result.json --matrices work/research_cycle8/context_design_audit/inventory_v1/nested_matrices.json --mapping work/research_cycle8/context_design_audit/exact_v1/mapping.json --out NEW_OUTPUT_DIRECTORY
```

Replay requires the exact pinned Python, NumPy and SciPy runtime. The result's scientific evidence can also be checked using standard-library rational arithmetic alone, without rerunning optimization: regenerate forced pairs from the stored intervals and matrix rows, verify each feasible beta and monotone completion, or verify each nonnegative normalized dual cancellation. No network or Git operation is needed for numerical verification.
