# Independent review: cPRIME exact transitions and anonymous-label assignments

**Verdict: the held numerical result and its stated finite-model proofs pass.** No mathematical or arithmetic mismatch was found. The principal qualification is interpretive precision at the twofold boundary and the distinction between recorded-point assignment robustness and continuous-envelope inference. A separate execution-guard recommendation appears below.

This is a separately authored review with inherited project context, not a blinded external review. No production numerical function was imported, no new source dataset or unselected phenotype was opened, and no endpoint, group, threshold or model was changed.

## Complete reconstruction

The independent standard-library checker made 6,350 assertions. It reconstructed all 60 group/cell profiles and all 720 assignments, checked all 2,160 within-assignment dispersion values exactly, and compared every profile against all 240 original widening-grid records. It also verified all 20 original exon common-order statuses. The maximum rational discrepancy is zero.

For the 48 two-marker cases, an augmenting-path bipartite maximum-matching implementation independently computes the lower and upper count bounds. Its results agree with full assignment enumeration in every case. For the 12 three-marker cases, all 36 assignments are retained and independently reconstructed. The matching implementation also agrees with exhaustive bijections on all 512 binary three-by-three graphs.

| Recorded-point statement at the fixed twofold threshold | Verified group/cell cases |
|---|---:|
| All three matched labels reach twofold under every admissible assignment | 8 of 60 |
| At least one matched label reaches twofold under every admissible assignment | 40 of 60 |
| The number reaching twofold changes with assignment | 23 of 60 |
| No assignment produces a matched label reaching twofold | 13 of 60 |

These categories overlap; their counts must not be summed. The 60 cases and 720 candidate mappings are not independent biological experiments. The summaries concern the fixed recorded point sets and the specified weak local assignment family.

## Proof challenge and verdicts

1. **Finite-marker transitions: correct on the stated positive, finite, closed interval domain.** Let A=max(L)/min(U) and B=max over distinct markers of U_i/L_j. There are at least two markers, and B>=1 follows because one of U_i/L_j and U_j/L_i is at least one. Thus the upper-dispersion formula is exactly B*t, without a missing clamp at one. Equality becomes possible at t=max(1,A), including the boundary. The lower-at-least-two and upper-at-most-two intervals have the inclusive endpoint conventions stated in the derivation. The inspected source has no incomplete interval records. The result must not be generalized to unknown or open endpoints without the separately stated extended conventions.

2. **Common-order boundary: correct.** For scaled intervals [L,tU], strict separation forces i before j precisely when t<L_j/U_i in at least one context. A topological order permits the prefix maximum of lower endpoints as a positive monotone completion. Any upper-bound violation would force a backward edge. Therefore acyclicity is necessary and sufficient for an arbitrary common scalar order with context-specific nondecreasing responses. At equality an edge is absent, which is essential to attain the boundary.

   The reviewer used max-min Floyd closure, not the author's sorted-weight scan. A closed walk's bottleneck cannot exceed the best simple cycle it contains, so the largest diagonal closure value equals the maximum directed-cycle bottleneck. Taking its maximum with one reproduces all five boundaries. Forty-eight invented systems also agree with exhaustive total-order optimization. The independently verified source cycles establish infeasibility at every t below the reported boundary; rational response witnesses establish feasibility at the boundary itself.

3. **Assignment gauge and counts: correct.** Fixing the first marker's labeling removes a simultaneous permutation of all labels. It leaves exactly 6 or 36 representatives for distinct row identities, including cases with coincident numerical values. Within-row RNA/DNA pairs are preserved. For two markers, good edges have dispersion strictly below two and bad edges have dispersion at least two. Completing a maximum good matching to a full bijection proves minimum bad count = 3 minus maximum good matching size; maximizing a bad matching proves the upper count. Inclusive twofold ties are handled on the bad side.

4. **Non-implication: retained and necessary.** The invented triples (1,4,16) and (2,8,32) have three twofold comparisons under every recorded-row assignment, yet their continuous envelopes overlap and permit equal responses. Assignment robustness does not increase calibrated confidence in a biological effect. Nor does dispersion establish a consistent signed marker effect. The author correctly separates this from the original envelope inference.

## Exact common-order boundaries

| Exon | Feasible from gamma | Certificate |
|---|---|---|
| 17 | [1.116006825003, 1.116006825004) contains the exact boundary | Two-edge cycle below; exact rational completion at equality |
| 25 | 1 exactly | Exact completion at the domain boundary |
| 27 | [1.051544002416, 1.051544002417) contains the exact boundary | Two-edge cycle below; exact rational completion at equality |
| 39 | 1 exactly | Exact completion at the domain boundary |
| 42 | 1 exactly | Exact completion at the domain boundary |

The exact t=gamma-squared values, every closure diagonal, all cycle origins and boundary responses are preserved in the reconstruction. Square-root display brackets were checked by rational squaring, including their strict upper endpoints. These are exact transitions under a declared hypothetical widening family, not estimates of measurement error. They concern unrestricted shared order and cannot establish additive representability.

## Decimal precision at twofold

Five group/cell cases contain assignment dispersions within 1e-12 of two. Their exact offsets are about 1e-15 to 1.6e-14. The 1e-12 value is a display-sensitivity diagnostic only; no classification was changed. All affected IDs and exact offsets appear in `exact_v1/reconstruction.json`.

The clearest interval-boundary case is `chr17-31225100-A-T`, KBM7. Its B value is exactly

```
2 + 1/468840859893519
```

so the maximum-at-most-twofold claim fails by approximately 2.13e-15 on the stored decimal inputs. The transition interval is therefore empty on t>=1, as the author reports. Treating the displayed value as exactly two would change that statement. This is numerical representation sensitivity, not evidence that the physical response exceeds a sharp twofold threshold.

The eight all-three-under-every-assignment cases do not depend on these near-ties. Across those cases, the smallest exact excess above two is `562427653717511781797930407/14303056614334166591159647955`, approximately 0.0393222. This separation does not supply a biological error bound; it only distinguishes them from the machine-precision boundary cases.

## Execution qualification

The held `analysis.py` uses Python assertions for the source hash and several certificate checks but does not reject `-O` or `-OO`. Its recorded normal execution is numerically verified here. Before offering the author script as a guarded archive entry point, add a reviewed explicit optimization refusal or replace those guards with explicit exceptions. This is not a defect in the inspected numerical results or a reason to recompute them. The independent replay checker already rejects optimized execution and checks both exact input hashes with explicit errors.

No broader novelty judgment is made: the protocol itself accurately treats interval orders, matching, assignment enumeration and necessary/possible reasoning as established methods. The useful contribution here is the complete source-specific transition and missing-label accounting, conditional on the deposited ratios and their measurement contract.

## Immutable review point and replay

- B1 source result: `eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf`.
- Reviewed follow-up result: `851bbf75b66dd0186d0b99d7536409773de1afdae5c3ea83c57563655d55308e`.
- Reviewed analysis source: `45f75f86ef97325ae17bcabd306504dc57f5d885e81d3244a71eb43a5870619f`.
- Reviewed protocol: `a3acdaad3beedba8793352fe402a587850c66a92bfb6c33b7ee42ca539c1decd`.
- Reviewed derivation: `fda4f419b8f6e3b906a58bba93bad31a94dca7130f67356ea57ef63fdc8f00d6`.

Run `verify.py` with explicit `--b1`, `--result`, `--protocol`, `--derivation`, `--analysis` and a fresh `--out` directory. It uses only the standard library, performs no network or Git operation, and does not import author code. The result and source hashes are checked before parsing and again after calculation. Run `decimal_boundary_check.py --reconstruction exact_v1/reconstruction.json --out NEW_DIRECTORY` for the separate precision receipt. Original author artifacts remain unchanged.
