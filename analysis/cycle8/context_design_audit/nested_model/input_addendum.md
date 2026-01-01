# D1 input and execution addendum

This addendum supplements the held design protocol, SHA-256 `53d9cd73e405a94741b9e25ef00bfe81f9cb4bc5e2a34776d0e6fa29f0b5f9c8`, without modifying it. Root authorized implementation and one execution after the implementation freeze and synthetic checks. Source outcomes were previously exposed under B1; this remains exploratory reuse.

## Exact inputs

- Existing B1 result: `work/research_cycle8/cprime_results/primary_v1/result.json`, SHA-256 `eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf`.
- Fixed matrix: `inventory_v1/nested_matrices.json`, SHA-256 `665c4d0586310d09868b07f01163875e82708c3441201ab7757d7c897b3ab874`, relative to the parent design-audit directory.
- Fixed state map: `exact_v1/mapping.json`, its exact SHA-256 recorded in `freeze.json` before interval access.

The source must have status `COMPLETE_DESCRIPTIVE_AUDIT`. Decode only the following fields: `state_cells` records' state, cell, usable_rows, complete, lower, upper and source_rows; source_rows records' source_row, state, cell and exon; sensitivities records' gamma and exon_shared_order_results' exon and status. Other fields are lexically skipped without JSON value decoding, including RNA, DNA, ratio, median, marker-group summaries and prior certificate details. No workbook, deposited counts or other source result is read.

Require all 132 state/cell interval records and all 396 source-row identities, exactly three per state/cell. Check canonical IDs and exon against the fixed map, source-row uniqueness and membership, exact matrix row order, model names and the existing gamma grid. Complete records require three usable rows and closed finite positive intervals. Incomplete records require the existing `(0,infinity)` representation, encoded as lower `"0"`, upper null. Exact finite bounds are rational strings; floating bounds are refused. No outcome-based exclusions or re-estimation occur.

## Implementation and numerical contract

`runner.py` uses the pinned matrices without feature selection. It generates only within-exon, within-cell forced comparisons. M0 shares coefficients across cells; M1 and M2 use the already frozen block matrices. Separate response functions are completed and checked for each context, including tied latent scores and the strictly positive domain of unknown intervals.

SciPy HiGHS supplies candidate feasibility or infeasibility information. Candidates alone are never results. Feasible coefficients are converted to rational numbers with denominator cap 1,000,000, then rescaled using their exact minimum positive forced gap. If the conversion cannot produce strictly positive exact gaps, the result is unresolved and execution stops. The cap is a numerical reconstruction setting, not an effect-size threshold or a search over scientific hypotheses.

For an infeasible candidate, a second linear program seeks a convex combination of constraint rows equal to zero. Support entries above 1e-9 are reconstructed with exact rational elimination. Free support variables are set to zero, and every weight, cancellation and normalization is then checked exactly. A failed reconstruction stops as unresolved; it does not authorize changing the model, threshold or source rows. No support threshold makes a scientific decision without exact certification.

Every feasible record saves coefficients, latent values, completed positive responses, all forced pairs and all exact margins. Every infeasible record saves all forced pairs and its normalized nonnegative rational dual. Independent verification should regenerate the rows and intervals directly from pinned inputs. The complete cases are 60 joint tests plus 80 individual M1/M2 context tests. Gamma nesting, model nesting, block conjunctions, saturation expectations and consistency with the existing weaker common-order test are mandatory checks. Any structural contradiction stops as a bug. A checkpoint preserves completed and unresolved records; the final result is written only when every check passes.

The runtime is CPython 3.12.14, NumPy 2.5.4 and SciPy 1.18.1. Versions must match before interval decoding. All implementation/document/source hashes and a successful synthetic receipt must be sealed in `freeze.json`; `authorization.json` must name that exact freeze and permit exploratory fitting. The runner refuses optimized Python, a pre-existing output directory, source drift, changed code or runtime mismatch. No Git or network operation is performed.

## Synthetic evidence before observed execution

`test_runner.py` passed 566 assertions. A separately implemented finite gap-cone oracle checks all 256 four-state interval combinations from closed points, overlaps and unknowns, and all 24 strict orders. Further checks cover exact witness/dual tampering, row and column permutations, reference recoding, shared versus separate cell order, tied and unknown completion, malformed matrix/interval inputs, selective-field poison guards, all 140 records on an invented full-size input, duplicate identities, source hashes, optimization and overwrite refusal. This is an independent algorithmic comparison within the same authoring lane, not a separate investigator or independent biological validation.

Observed execution is a single fixed run after this addendum and implementation are sealed. Root will independently inspect the proof and reconstruct numerical certificates before admitting any result into the research record. Failed or unresolved results remain failures. No positive outcome is required.
