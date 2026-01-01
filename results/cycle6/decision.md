# FAS: a measured ordering constraint, without a certified barrier

The deposited included-RNA enrichment measurements exclude the declared codon-additive latent score followed by an arbitrary common nondecreasing response **within the prescribed empirical envelopes**. All 72 variants in the native-protein class of the published design have finite envelopes. An exact five-comparison cancellation proves the exclusion. Its narrowest relative separation is 1.36%, which is a margin of this particular witness, not the robustness radius of the complete system.

No accessibility restriction is certified across the primary envelopes for any of the 2,556 endpoint pairs. The corrected point summary has 73 barrier pairs, and individual normalization scenarios have 62 to 307. This difference is a substantive uncertainty limitation. It neither proves the unknown biological landscape connected nor establishes that any particular point barrier requires interactions.

The evidence is credible as a finite computational statement: separate implementations reconstruct every selected count, normalization ratio and envelope, and independently verify all exact certificates and path bounds. Its biological interpretation remains limited by genotype-dependent RNA abundance, stability and recovery. Those effects can mimic an ordering that a common monotone codon-additive model cannot represent. No independent-study validation or calibrated biological uncertainty is supplied by the repeated transfections.

## Three distinct questions

| Question | Verified answer | Quantitative evidence | Limitation |
|---|---|---|---|
| Does the declared model fit inside the envelopes? | No. Exact conditional exclusion. | 74 forced orders; five rows with equal weight 1/5 cancel exactly. | An exclusion of a representation of this assay endpoint, not direct molecular interactions. |
| Is a restriction certified in the full declared graph? | No pair has a nonempty conservative restriction interval. | 72 vertices, 204 edges, 2,556 pairs, zero unknown states, zero certified pairs. | Coordinatewise envelopes are conservative sensitivity sets. Unknown biology need not realize every combination. |
| What connects the two results? | No barrier-specific interaction attribution follows. | Point and single-scenario barriers coexist with ordinal incompatibility. | Co-occurrence does not show that the particular barrier requires departure from the model. |

## Complete result and source reconstruction

The source is Baeza-Centurion et al., *Cell* (2019), DOI [10.1016/j.cell.2018.12.010](https://doi.org/10.1016/j.cell.2018.12.010). The twelve deposited files are pinned to `lehner-lab/Scaling_Law` revision `b87064c39a553b0759d86c2b428813c58dc17ae2`. All original Git blobs were checked. Nine outputs are separate transfections of one library; three DNA inputs are technical PCR preparations. The 27 cross-combinations share those preparations.

The native included-protein annotation uses the full NM_000043.6 coding frame and NP_000034.1. It does not authenticate translation of the experimental minigene or preservation of isoform mixture. The 72-state graph is complete within the declared 3,072-member design, not within the ambient synonymous sequence space.

All 72 selected variants have usable records in every source file. Four output files each omit one design record outside this class. No unknown genotype was removed to improve connectivity. Whole-file and TR1-aligned normalization totals agree. Every primary and secondary endpoint-pair comparison is retained in the full result.

The original parser stopped at `Unexpected count header` before parsing outcomes. An independently reviewed amendment restores the expected column names in memory and preserves every original row. Scientific freeze `7d8faf59` and format amendment `9ed05ad8` were committed before successful evaluation. Source examples and the format-diagnostic reference rows were already exposed and are logged. The public R import would consume the native-reference row if applied literally to the headerless deposited files. Consequently, the secondary point is a **corrected headerless source-formula reconstruction**, not authenticated reproduction of the historical point estimate. This observation does not establish an error in the authors' actual working files or published conclusions.

## Exact certificate

The fixed zero-based indices satisfy the five strict empirical orders
`16 < 46`, `35 < 46`, `39 < 46`, `42 < 51`, and `70 < 13`.
For the eight-column categorical design, their difference vectors sum to zero:

`(X46-X16) + (X46-X35) + (X46-X39) + (X51-X42) + (X13-X70) = 0`.

A common nondecreasing response would require every corresponding latent difference to be strictly positive. Their zero sum makes this impossible. No optimizer status alone is used as proof. Exact rational bounds, weights and feature differences are provided in `FAS_Exact_Order_Certificate.csv`.

The five separation ratios are approximately 1.12909, 1.22291, 1.10198, 1.01365 and 1.03673. The smallest absolute gap is 0.0075674 relative-enrichment units. No optimality or biological calibration is claimed for these witness margins. The source already reports synonymous context effects and specific interactions at overlapping positions. This analysis adds a transparent full-domain, conditional model certificate, not their first biological discovery or a new foundational theory of order representation.

## Connections, restrictions and uncertainty

For every endpoint pair, the tables give the lower-field and upper-field maximum path bottlenecks, `C_L` and `C_U`, and the minimum endpoint lower score `E`. A connection is guaranteed for thresholds at most `C_L`. A restriction requires a nonempty interval `(C_U, E]`. All other threshold cases remain unresolved or lack guaranteed acceptable endpoints. All 2,556 primary restriction intervals are empty.

The point has 73 nonempty intervals. Every one of the 27 dependent normalization fields has an exact ordinal contradiction and between 62 and 307 barrier pairs. These secondary findings are preserved. They are not independent confirmations, confidence levels, or a replacement primary endpoint. The rectangular sensitivity box can combine extrema never observed together. Testing observed joint scenarios is a separately labeled exploratory follow-up, not a revision of the zero primary count.

## Independent checks

- Independent source reader: 1,944 exact ratios, 648 per-output relative scores, 72 means and 144 envelope endpoints agree with zero discrepancy. The total selected-field comparison count, including identities and count entries, is 3,816.
- Independent mathematical implementation: all 30 rational certificates and 74,124 primary/secondary pair records agree, using descending-threshold graph searches rather than the production closure algorithm.
- The replay works in a relocated directory without Git metadata or network access. Original verifiers, the format failure and review, frozen code, complete source counts and compressed full result are retained.
- The figure displays all 72 genotypes, all 27 scenarios and the exact returned certificate. Exported figure boundaries and final-size typography were independently inspected.

## Single most informative biological test

The key missing record is a separately prepared, sequence-linked inclusion-and-skipping measurement for all eight certificate variants, with matched reference controls, preparation identifiers and a defensible differential recovery/error calibration. A pooled skipped transcript loses the exon bases that identify its genotype, so a linked barcode or individually identified clone is necessary to assign that denominator. A further enrichment-only measurement does not resolve the same ambiguity.

The earlier FAS source is currently ineligible for this exact discrimination test: only five of eight certificate vertices are in its analyzed single/double-mutant design, and its direct clone validations concern single mutants. None of the seven non-reference certificate variants is a single mutant. These are design memberships, not observed coverage. No earlier FAS outcomes or presence table were opened, and no favorable subset was substituted.

**Continuation decision:** close this primary test without retuning. Proceed only with two explicitly exploratory computations justified by the outcome: the full model's perturbation threshold and finite-joint-scenario accessibility. Preserve the orthogonal measurement gap, retired Cambray validation arms and all failed predictions. Existing published releases remain unchanged. A new release has not been completed for this cycle.

The full uncompressed result has SHA-256 `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466`.
