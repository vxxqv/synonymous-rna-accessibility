# Independent FAS mathematical result audit

Verdict: PASS for the saved finite-domain rational certificates, ratio arithmetic and graph bounds. No disagreement was found. This is not an independent source-count extraction, a biological replication, or an inferential calibration.

The explicitly authorized input `full_result_v1.json` has SHA256 `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466`. The independent verifier was fixed and synthetically tested before its observed execution. It imports no production implementation. All nine original freeze-listed files still match their stored hashes, and the result's freeze hash matches the unchanged freeze file. The later header-format amendment and historical access ordering are not independently authenticated by these arithmetic checks.

## Exact result

- 72 declared vertices, 204 independently reconstructed nucleotide edges; 72 finite envelopes, zero unknown envelopes.
- 1,944 scenario ratios, 648 author-relative output values and all 72 point means agree exactly with reconstruction from the saved count vectors and supplied TR1-aligned totals.
- All 30 rational certificates validate: primary codon, equivalent nucleotide encoding, 27 normalization scenarios and the author-relative point summary. All 30 report infeasibility and have valid exact nonnegative cancellation witnesses.
- The primary system has 74 forced strict orders. Its supplied certificate uses five orders and eight distinct vertices.
- Independent descending-threshold BFS components agree with all 74,124 saved pair records: 2,556 primary records and 28 times 2,556 secondary records.
- Primary conservative obstruction count is **zero**. The individual normalization scenarios have between 62 and 307 obstructed pairs; the descriptive point field has 73. These dependent scenarios are not independent confirmations, and their barriers do not overturn the primary zero count.

## Compact primary certificate

Vertex indices are zero-based in the frozen result. The following five observed envelope separations are strict:

| Lower vertex i | Higher vertex j | U_i < L_j gap, approximate |
|---:|---:|---:|
| 16 | 46 | 0.11432955561800273 |
| 35 | 46 | 0.18227923721130310 |
| 39 | 46 | 0.09254565969008250 |
| 42 | 51 | 0.00756741397787093 |
| 70 | 13 | 0.02083279918361126 |

Exact rational bounds, gaps, weights, codon states and differences are in `primary_certificate.csv` and `observed_v1/verification.json`. The rounded values above are explanatory, not used in acceptance tests.

Reconstructing the eight-column codon design directly from the frozen sequences gives the exact identity

`(X46-X16) + (X46-X35) + (X46-X39) + (X51-X42) + (X13-X70) = 0`.

Any common nondecreasing response consistent with the five disjoint envelope pairs would require each corresponding latent difference to be strictly positive. Their sum must instead be zero. Equivalently, after the purely mathematical unit-margin rescaling, averaging these five inequalities with weights 1/5 gives the contradiction 0 >= 1. Thus no arbitrary codon main effects followed by a common nondecreasing function can place all 72 fitted values inside these empirical envelopes. No numerical optimizer verdict is needed once these exact identities and strict gaps are checked.

This is an ordinal descriptive constraint, not a measured interaction magnitude or a p-value. It rejects that common monotone codon-additive representation of the declared source endpoint within the stated envelopes. It does not prove a molecular interaction, a splicing mechanism, robustness to unobserved measurement error, or biological function preservation. Variant-dependent recovery bias remains an alternative. The unit latent gap is not a biological effect-size threshold. The equivalent nucleotide certificate is separately checked, and its different support is not a separate biological finding.

## Graph interpretation

The independent algorithm activates all vertices at each distinct exact score, including ties, then uses BFS to form threshold components. The first level at which each pair coalesces is its widest-path capacity. This differs from the production max-min Floyd closure. Every selected witness includes an explicit path at the capacity and the endpoint's component at the strictly greater level, or records endpoint blockage at that level.

For lower and upper fields, each reported robust threshold interval is exactly `(C_upper, min(L_s,L_t)]`; it is nonempty only when the left endpoint is strictly smaller than the right. Equality gives an empty interval, although acceptance at a given threshold remains weak. Positive infinity for an unknown upper value leaves that vertex optimistically available. Disconnected base endpoints are handled by negative-infinite capacity in the synthetic verification; the observed 72-vertex graph is connected.

The primary zero count means that this certificate does not force any endpoint-relative obstruction across all fields inside the preparation envelope. It does not establish connection for every admissible field, absence of particular-scenario bottlenecks, or absence of barriers elsewhere in the unmeasured synonymous graph. Null rejection and barriers are separate claims; neither attributes a particular bottleneck to causal interactions.

## Verification coverage and limits

Before observed execution, all undirected graphs on at most four vertices with scores from {-1/2, 0, 1, infinity} were compared with an independently implemented exhaustive simple-path oracle: 16,932 graph/score cases and 167,012 capacities. Tests also cover exact feasible certificates, constant responses and latent ties, unknown upper bounds, touching intervals, altered dual weights/differences/orders, and rejection of -O/-OO before input access. Existing output is rejected before input parsing. These finite tests support implementation correctness, not a proof for arbitrary graph size; the threshold-component argument supplies the general characterization.

No raw count file, new source, model, threshold or subset was opened or selected for this audit. The full 3072-row normalization totals are supplied metadata, not independently recomputed here. Missingness rules are tested synthetically but not empirically exercised, since all observed envelopes are finite. The saved 72-row counts, their ratio identities and mathematical consequences are checked. Shared-library and shared-denominator dependence is preserved; empirical envelopes remain descriptive rather than confidence bands. No scientific threshold was adjusted.

## Replay

Run from the workspace root using Python 3.12 or compatible standard-library Python:

```
python work/research_cycle6/fas_result_math_audit/verify_result.py --result work/research_cycle6/fas_results/full_result_v1.json --result-sha256 3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466 --contract outputs/synonymous-rna-accessibility/analysis/cycle6/fas_pilot/source_contract.json --contract-sha256 f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38 --out NEW_AUDIT_DIRECTORY
```

The output directory must not exist. A new synthetic receipt can be produced by copying the two synthetic/verifier source files into a fresh directory and running `test_synthetic.py`; the original saved receipt is never overwritten. Full hashes are in `receipt.json` and `synthetic_receipt.json`.
