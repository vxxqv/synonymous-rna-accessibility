# Independent FAS count and normalization audit

**Verdict: PASS for the count, normalization and eligibility calculations at the exact saved revision.** The independent implementation uses only the Python standard library and exact rational arithmetic. It does not import the production analyzer or reuse its intermediate numerical values. Observed input access began only after the format amendment was committed and the primary run was authorized.

## Exact scope and evidence

The audit verified amendment commit `9ed05ad8ab75c596e87a5ce7b63eb1bb38e71ba5`, the scientific and format freeze manifests, every frozen implementation hash, the 72-sequence source contract, and all 12 original count-file Git blobs before decoding count rows. It retained the first native-reference row, checked every source row for valid identity, integer counts, uniqueness and the ancestral mutation annotation, and rechecked source and result hashes after the calculation.

All nine output files and all three technical-input files were read. Eight files contain 3,072 observed rows; BR4, BR5, BR7 and BR9 contain 3,071. Each of the 72 selected genotypes has all required counts, so these four missing rows do not cause selected-domain exclusions. All 72 empirical envelopes are finite, with zero unknown envelopes. This is complete coverage of the fixed selected class, not evidence that the entire library has complete measurement coverage.

| Independently reconstructed quantity | Number checked | Result |
| --- | ---: | --- |
| Selected integer count entries | 864 | Exact agreement |
| Output/input relative enrichment ratios | 1,944 | Exact rational agreement |
| Empirical lower and upper endpoints | 144 | Exact rational agreement |
| Corrected per-output source-formula points | 648 | Exact rational agreement |
| Corrected mean source-formula points | 72 | Exact rational agreement |
| Scenario unknown-vertex counts | 28 | Exact agreement |
| TR1-aligned totals, whole-library totals and normalization metadata | Complete saved object | Exact agreement |

The 3,816 selected-field scalar comparisons include 144 identity/index comparisons in addition to the numerical fields above. The maximum exact numerical discrepancy is zero. No failed or undefined selected cases were suppressed.

## Separate reconstruction

For every selected genotype, the reference calculation divides each output/input count enrichment by the corresponding native-reference enrichment. Its 27 entries retain the fixed order BR1/TR1 through BR9/TR3. A missing required count or a required zero denominator produces an unknown value; a zero focal output count remains an observed zero. Complete envelopes are the minimum and maximum of the 27 exact ratios.

The secondary calculation explicitly restricts every library total to the TR1 genotype universe. It calculates the median of available normalized technical-input frequencies, then the normalized output enrichment relative to the native reference for each of the nine outputs. The reported scalar is the mean of those nine relative enrichments. This agrees exactly with the saved result. The point lies within its complete empirical envelope for all 72 genotypes.

The secondary quantity is a **corrected headerless source-formula reconstruction**. It is not an authenticated reproduction of the historical author result: the deposited first native-reference row is retained here, whereas applying the original header-setting literally would drop that row. The saved legacy scenario label does not change this interpretation.

## Synthetic safeguards

The separate synthetic receipt records ten checks. Known-answer values distinguish whole-library totals from TR1-aligned totals, verify all 27 ratios and the exact median calculation, and cover a missing technical input, zero technical input, zero focal output, missing output, zero native-reference output and missing TR1 identity. In particular, a missing denominator produces an unknown envelope and a zero focal output produces a valid zero. These checks were run before the observed comparison.

## Limits and exposure

The 27 scenarios reuse nine outputs and three technical inputs. They are dependent preparation combinations, not 27 independent biological replicates, and their empirical envelope is not a calibrated confidence interval. This audit does not independently recalculate the ordinal feasibility decision, the 204-edge graph, or the barrier certificates; those require the separate geometry audit. Exact agreement checks the specified reconstruction, not an unverified source measurement interpretation.

The verifier decoded all 12 pinned count files to obtain the TR1-aligned totals and checked the entire fixed 72-genotype class. Only compact audit metadata are reported here. No new model, subset, threshold, exclusion or outcome-dependent selection was introduced. This review is separately authored with inherited project context, not blinded external review.

## Reproduction and fixed hashes

Run from the project workspace with the existing Python runtime:

```text
work/env/Scripts/python.exe work/research_cycle6/fas_result_numerical_audit/verify_counts.py --mode observed --authorized-observed --repo outputs/synonymous-rna-accessibility --input-directory work/research_cycle6/fas_inputs --result work/research_cycle6/fas_results/full_result_v1.json --out NEW_RECEIPT.json
```

The output must be a new path. The explicit observed flag records the caller's authorization; it does not replace authorization. This verifier uses local Git to check the committed amendment and makes no network request.

| Input or evidence | SHA-256 |
| --- | --- |
| `verify_counts.py` | `f5a932166d3cd8ddefc47020e190610b0b845aea9a72a1d6a3531f9693ffd840` |
| `synthetic_receipt.json` | `47d7d57f65d50f774d5ac6b21a6e948e2b3e2d1236b06afe124c1e66002af36c` |
| `observed_receipt.json` | `b241c91c1b639b4753eacbc119d32d99aadb1e7ad57b1b9ac92e2a258accad77` |
| `full_result_v1.json` | `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466` |
| `source_contract.json` | `f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38` |
| Scientific `freeze.json` | `7bddb75dcdb4b401a120dcb518bd34fcfb5f6a801a492ee1a59f1fced90be3eb` |
| `format_freeze.json` | `d3b69ed44015a712d4dbe436cdb7399c009e521f3444a6a5736c547a04094041` |

Per-file source SHA-256 values, Git blob identities and row counts are retained in `observed_receipt.json`. There are no outstanding numerical discrepancies within this audit's scope.
