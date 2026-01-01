# Independently authored source-first reference successor

This version is a post-DP-exposure source-contract correction, sealed before actual CAO, normalized assays, contrasts, predictions or performance were opened by this reviewer. It is not the original prospective reference. The original reference.py and preaccess_seal.json remain unchanged. Production v3 was reviewed before this correction; the successor is not a blinded independent rediscovery of its depth rule.

The sole mathematical change is replacement of the unique-depth assertion with depth supplied by the lexically first collapsed raw VAR_ID, within the reviewed fixed-sample/ASCII source contract. Exact per-allele median collapse, total alternate counts, WT subtraction, NC handling, logs, missingness and preparation summaries are retained. The comparator additionally authenticates selected supplier, ordered depth alternatives, per-allele depths, conflict flags and the source-rule label. Negative reconstructed WT stops immediately, records sample/position and exact selected-depth/total/WT details, and attempts no repair.

All 25 other common functions are byte-identical to the original reference, including every regression/isotonic, prediction, fold, scoring, sign and decision function. The only changed common functions are normalize, compare_normalization, read_tsv and run. New functions are provenance and guarded-comparison helpers. read_tsv now decodes only needed identity, DP/CAO and mask fields; unrequested CAF and RF-feature bytes remain undecoded. This is a narrower reader, not a new scientific filter.

Synthetic receipt: 288 original-helper comparisons across 24 invented cases/orders, lower-first versus maximum counterexample, negative-WT stop and five deliberate lineage falsifiers pass. All 2,500 predictions, 2,500 errors, 80 aggregate scores, 3,160 site errors, 40 coefficient/rank sets and 500 range-flag sets still match the earlier invented stage exactly. Complete downstream objects equal the original reference's objects. Twenty-eight metadata-only authorization/lineage tests pass. -O/-OO and overwrite are refused. No actual scientific numeric payload was decoded for these checks.

## Actual comparison requires separate authorization

```text
python reference_source_first.py --authorization EXPLICIT_AUTHORIZATION.json --out NEW_DIRECTORY
```

Authorization schema: COMT_independent_source_first_authorization_v2. Set status explicitly_authorized and reference_sha256 to the sealed file hash. inputs maps every role below to exactly path and sha256. The observed comparison requires all 21 roles:

| Role | Artifact |
|---|---|
| protocol | Original protocol_v2.md |
| integration_addendum | Unchanged integration_addendum.md |
| normalization_amendment | normalization_contract_amendment_v3.md |
| predecessor_reference | Unchanged original independent reference.py |
| predecessor_reference_seal | Unchanged original preaccess_seal.json |
| predecessor_code_manifest | Original production code_freeze_v1.json, SHA f63371... |
| successor_code_manifest | Production v3 code freeze |
| contract | The unchanged 125-target/53-site metadata contract |
| extract_manifest | Original sealed extraction manifest |
| rf_manifest | Original sealed RF manifest |
| normalization_manifest | Successor normalization seal |
| prediction_manifest | Successor prediction seal |
| evaluation_manifest | Successor evaluation seal |
| rows_1, rows_2 | Explicit extracted biological TSVs, not source simulation inputs |
| row_status_1, row_status_2 | Explicit frozen RF mask TSVs |
| normalization | Successor normalization.json |
| contrasts | Successor contrasts.json |
| predictions | Successor predictions.json |
| evaluation | Successor evaluation.json |

All hashes and stage/predecessor associations pass before numerical bodies are decoded. The stage seals must bind the exact selected payloads and six stage links. The normalization amendment, original reference/seal and predecessor production freeze are hard-pinned. An explicitly marked synthetic_fixture_only authorization may omit the five raw roles for wrapper tests; an observed authorization cannot silently downgrade full reconstruction.

Count and model fractions must agree exactly; unchanged relative/absolute tolerance 2e-12 applies only to log/preparation arithmetic identity. Labels, signs, missingness, coverage and scientific thresholds are exact. Independent counts and summaries are derived before production answers are opened; the verified stored contrasts then supply the specified stored numeric inputs to exact model arithmetic. A substantive mismatch writes a failure receipt and stops without retuning. No Git/network/install or original-workspace path is needed at runtime. Destinations must be fresh.

Limits: extraction and frozen RF masks are conditioned inputs, not independently trained/classified here. The checker verifies the selected payload/hash lineage, not every unselected stage file, biological sample identity, historical access chronology, physical WT meaning or calibrated uncertainty. It does not reproduce model condition numbers. Agreement cannot establish protein amount, mechanistic compensation, an independent load effect or assay equivalence.
