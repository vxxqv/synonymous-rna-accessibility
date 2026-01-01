# Operational successor adapter preparation

Status: ready for independent read-only review. No actual grouping, fitting or new outcome access has occurred. Root must authorize and freeze the new protocol and all source evidence before these adapters can be used. This directory is separate from the untouched v1 repository and original work files. Exact original bytes are retained in `originals/`.

The new source-gate schema is `CNOT3_operational_source_gate_v2`, with status `operational_cohort_verified`. It must state `operational_source_gate_passed=true`, `original_frozen_gate_status="held_for_source_diagnosis"`, `historical_published_cohort_reproduced=false`, and exactly `unresolved_source_discrepancies=["published_cohort_size_mismatch"]`. It must not set `reproduction_gate_passed=true`.

All11 operational checks must be present, true and have exactly these names:

```text
count_5480
distinct_stable_ids
positive_finite_paired_half_lives
positive_finite_SE_IP
STAT_IP_arithmetic
no_negative_SE_among_IP_candidates
no_top1000_boundary_tie
median_log2_change_top1000_greater_rest
ten_significant_digit_ratio_compatibility
all_49_original_ratio_discrepancy_genes_retained
independent_row_join_verified
```

The gate carries `cohort_gene_count=5480`, the sorted `cohort_GeneIDs`, `cohort_rows_sha256`, `prediction_protocol_sha256`, `original_frozen_gate_sha256`, `independent_source_audit_sha256`, and sorted `original_ratio_discrepancy_GeneIDs`. The original49 IDs are extracted again from the supplied original gate and compared exactly with this list. All must belong to the unchanged5480 cohort. The new gate's original-cohort IDs and row hash must equal those in the preserved failed gate. The count defect and original ratio checks are not deleted or declared successful.

`pipeline.py` adds required input roles and CLI options `original_gate`/`--original-gate` and `source_audit`/`--source-audit`. Their actual file bytes are hashed with all other inputs. Root authorization must pin every input hash and the two code hashes, use `operation="CNOT3_operational_cohort_group_partition_v2"`, and set `real_execution_authorized=true` and `operational_source_gate_passed=true`. Unverified hash strings in a gate alone are insufficient. The original/source-audit hashes must match these supplied and authorized files before rows, reference selection or grouping proceed. The root-authenticated gate summarizes the separate source audit; this adapter does not redo that audit or infer an actual generating procedure from rounding compatibility.

The source row schema is unchanged. Reference hashes, stable-ID rules, eligibility functions, grouping algorithms, caps and hash salts are unchanged. The sole `cohort.py` change is the explicit design guard from `reproduction_count_5481` to `operational_count_5480`. The historical eligibility-report field `input_reproduction_genes` still refers to the number of rows received from the existing reconstruction, not a claim that it passed. This one-line guard change accounts for the one-byte file size difference from the original. No text-encoding conversion was substituted for a numerical change.

The downstream prediction design receipt uses `status="root_authorized_operational_prediction"`, `analysis_status="post_exposure_exploratory"`, `operational_source_gate_passed=true`, the false historical reproduction flag and held original status above. It requires SHA256 fields for the operational gate, original gate and independent source audit, alongside the unchanged reference/grouping success fields. The coverage denominator is named `operational_slam_overlap_count`; its mathematical definition and80% threshold are unchanged. The prediction freeze must have the same new authorization status. Evaluation accepts `root_verified_operational_prediction_commit`, with the same prediction/target/code hash checks as before. Source provenance is retained inside the frozen design receipt. No holding-out language erases the prior source reproduction exposure.

`VERSION` is now `cnot3-operational-exploratory-v2`. Every numerical prediction, residualization, tuning, transformation, fitting, scoring and influence function body is AST-identical to the reviewed original. All model constants and grouping/resource thresholds are unchanged. The accepted function-level AST differences are only: `cohort.py:design_guards`; `pipeline.py:execute,main`, plus new guard `verify_operational_gate`; and `predict.py:check_design,run_predict,run_evaluate`. The latter differences concern authorization/schema/status only. `predict` and `evaluate` themselves are unchanged. `exact_changes.diff` exposes every textual modification, and `ast_evidence.json` records hashes for every unchanged function/class.

Eight synthetic guard tests pass. They reject old gate schemas, false pass labels, missing/false/unknown checks, changed input/evidence hashes,5481 counts, changed IDs, omission or substitution of any of the49 retained rows, changed historical failure evidence, and prediction receipts that hide exploratory status. A modified synthetic source-audit file is rejected before references or grouping can run. Hash salts and thresholds are checked explicitly. Tests do not call grouping or fitting. The previous numerical synthetic tests need not be rerun merely to establish unchanged numerical ASTs; independent review can separately verify the adapters and existing numerical evidence.

Reproduce the bounded preparation audit using `work/env/Scripts/python.exe work/decay_cycle4/operational_successor/implementation/audit_preparation.py`. It produces a diff, AST record, guard-test log and preparation receipt only. Actual scientific phases must still use the separately preserved runtime/script guard for CPython3.12.14, NumPy2.5.4 and Biopython1.88. No actual phase authorization is included here.
