# Independent observed robustness result review

**PASS.** All 26 reported radii are supported by an exact boundary construction and an exact contradiction persisting at every smaller admissible parameter. All 72 observation-deletion certificates, all 12 preparation-omission sets and all 220 stored search-event certificates pass. Maximum exact discrepancy is zero. The verifier uses standard-library rational arithmetic and imports no production radius routine or optimizer.

## Verified full-domain result

The fixed native-reference identity remains exactly one. The sharp additive expansion is

`epsilon = 519286149121 / 31050521652448`, approximately `0.01672391` relative-enrichment units.

The sharp squared multiplicative parameter is

`t = gamma^2 = 9128860950384 / 8609574801263`.

Its exact decimal enclosure is

`1.029715960236 <= gamma <= 1.029715960237`.

These values agree with the saved result. They describe the minimum expansion of the declared empirical envelope admitting at least one common monotone codon-additive representation. They are not calibrated measurement errors, confidence bounds, molecular interaction strengths or evidence that every permitted perturbation repairs the model. A claim of biological robustness would require independently justified error or recovery bounds smaller than these thresholds.

## What was independently reconstructed

The verifier checked the pinned successor freeze, all its input source hashes, the compressed and decompressed primary result, the source sequence contract and the exact saved successor result. It reconstructed the eight-column codon design from the fixed sequence strings, checked the native identity at vertex 46, and recomputed every primary and preparation-deletion envelope from all required source scenario ratios. Source bytes were checked again after the audit.

For each of the two full-domain radii and 24 preparation-deletion radii, it independently enumerated every strict-order expiry from the input intervals and binary anchor mask. This verifies all 26 complete event inventories. It then checked each stored search certificate at its exact event, the boundary interval arrays, rational latent coefficients, common nondecreasing response including ties, and the immediately preceding exact dual. Each positively weighted dual order either persists permanently or expires no earlier than the proposed boundary. Positive weights sum to one and the feature differences cancel exactly. Thus no smaller admissible parameter is feasible, while the boundary is feasible. No reoptimization was performed.

For multiplicative results, the boundary response is in squared units. The independently checked nonnegative square-root transformation gives the physical response. Both rational gamma enclosure endpoints and their exact decimal strings were verified. Every native-reference boundary interval stayed [1,1].

All nine output deletions retained exactly the other 24 scenario columns; all three input deletions retained exactly the other 18 columns. Every resulting envelope was nested inside the full envelope and every radius respected the corresponding monotonicity. Six preparation omissions leave both radii unchanged; six increase them. The complete exact values are retained in the receipt and `preparation_radii.csv`. No preparation was omitted from reporting or chosen after inspecting its effect.

## Observation influence and the reference distinction

Exactly three formal single-row deletions make the remaining system feasible: vertices **13, 42 and 46**. The other 69 deletions retain an exact contradiction. The retained identity list and every local certificate index were checked independently.

Vertices 13 and 42 correspond to measured genotype observations in the fixed domain. Their deletion reveals dependence of this full-system exclusion on those individual constraints. This is a within-dataset influence result, not a claim of causal biological importance.

Vertex 46 is different: its response equals one by normalization. It is not an independently noisy measurement of the same kind. Deleting its observation row removes the fixed reference-identity constraint while all other rows still use their original reference normalization. The resulting feasibility must therefore be reported as a **structural relaxation of the normalization identity**, separately from the two measured-observation omissions. It must not be described as a third biological outlier, an independent noisy replicate, or evidence that recomputing the normalization would restore additivity.

## Limits

This is a retrospective exact sensitivity analysis of exposed primary data. The 27 joint preparation fields are dependent, and their coordinatewise box need not preserve that dependence. Removing preparation scenarios narrows the box by construction, so the observed nondecrease of radii is a mathematical consistency check, not an independent robustness experiment. The calculation neither attributes the primary exclusion to a molecular mechanism nor establishes an accessibility barrier. The original primary conclusions are unchanged.

The review is separately authored with inherited project context. It checks every reported finite-domain certificate and declared omission, but does not supply new measurements, independent biological validation or a model of unmeasured recovery bias. No result was retuned and no additional subset was selected.

## Reproduction and evidence

```text
python work/research_cycle7/robustness_result_review/verify_result.py --repo outputs/synonymous-rna-accessibility --result work/research_cycle7/results/robustness.json --out NEW_REVIEW_DIRECTORY
```

The output directory must be new. The verifier needs Python's standard library only and does not depend on Git or network access. It refuses optimized Python, drifted input hashes and existing outputs. The actual run took 1.79 seconds using CPython 3.12.14.

| Artifact | SHA-256 |
| --- | --- |
| Saved successor result | `861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa` |
| Exact primary result | `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466` |
| Successor freeze | `340624e795e600abbbaeb3daefa89b054ceca7e356855014d320df20b043a10f` |
| Independent verifier | `e971da5fab874f0843b58d2e85dc9fc6b28939254eb4ee3c19d73466f46bc46f` |
| Observed audit receipt | `e2be8de1c903efa19b4ad8e2f47161e16d76f2f5ae1541b1bb003f790efe4ecf` |

The declared freeze commit is `7eae1381c781c93e324f01d4e0d07692bceabea5`, whose local presence was checked. The numerical verifier authenticates frozen bytes rather than relying on commit dates. No arithmetic or certificate failure was found.
