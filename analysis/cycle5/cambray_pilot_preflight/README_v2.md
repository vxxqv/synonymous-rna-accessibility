# Shared-corner masking revision

The original preflight manifest and all eight files it pins remain byte-identical. `preflight_manifest_v1.json` is an exact copy of `preflight_manifest.json`. This revision adds a union mask and does not rewrite the original claims, split or receipts.

For each arm, the mask is the union of **all original** per-square hashed withheld IDs. It is formed before exclusions and is never recalculated after excluding a square. A square is prediction-eligible only if it contains exactly one member of that fixed union. Exclusions affect prediction only.

| Arm | Descriptive squares / diagonals | Union-mask IDs | Prediction squares | Related seed groups | Distinct prediction targets |
|---|---|---:|---:|---:|---:|
| Development |23 /46|22|22|16|21|
| Reserved validation |24 /48|24|24|16|24|

Development square45 in seed47 contains both `47_12111223_1` and `47_13111223_2` in the union mask, so it is excluded from prediction with reason `multiple_globally_masked_corners`. Its original target stays masked and both descriptive diagonals stay retained. Seed47 retains two prediction squares. No other square is excluded. Two eligible development squares share a target, which is additional outcome dependence; neither the22 development cases nor the24 reserved cases are independent experimental replicates.

All47 squares and94 diagonal pairs remain in the primary descriptive inventory. Both prediction arms satisfy the10-square/eight-group minimum at metadata level. Finite measurement coverage and the requirement for both true classes remain unexamined. The masking changes neither the model formulas nor chosen corners, groups, thresholds or success criteria.

`predict_from_visible_ids` accepts exactly the three unmasked values from one eligible square, in the fixed coordinate convention. It refuses extra table entries, a masked value, missing known corners or a square with multiple masked IDs before accessing values. It has no file reader or alternative table source. The future staged reader must enforce the hash-pinned union mask before decoding fields. Predictions must be saved and sealed before masked values are opened, including when a target is also a corner in another square. Target missingness is evaluated only after that seal, not used to obtain early target information.

Five new tests pass: overlapping-square leakage rejection; duplicate target handling; strict key validation with poison values; cross-arm masked identity refusal; and exact real metadata invariants. All original pinned hashes are rechecked. No new phenotype values are read. `prediction_mask_v2.json` includes every unit, exact masked sets, exclusions and counts by related seed group.

The root has also resolved the two earlier protocol details: author-mean comparison uses relative tolerance1e-8 and absolute tolerance1e-10 with all differences retained as rounding diagnostics; the classification screen compares seed-balanced error against each simple baseline and requires both classes plus at least10 prediction squares in eight groups. These choices are recorded here, not evaluated on outcomes.

Fresh-output-only reproduction commands:

```text
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_pilot_preflight/masking_v2.py
work/env/Scripts/python.exe work/joint_function_cycle5/cambray_pilot_preflight/test_masking_v2.py
```

These commands refuse to overwrite their own new result files. Existing v1 files are read and hash-verified, not edited. Separate review and protocol freeze remain prerequisites for real phenotype work.
