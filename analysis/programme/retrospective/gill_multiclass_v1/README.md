# Frozen Gill multiclass measurement protocol

This is a pre-outcome protocol, not a completed measurement result. It retains 2,176 complete connected local protein classes, with 1,052 classes reserved for the primary descriptive comparison. All belong to one closely related reporter system. The original reference-protein prediction failure remains unchanged.

`protocol.json` is the authoritative finalized protocol. `freeze_proposal.json` and `biological_requirements.md` retain the outcome-blind proposal. The primary agent's separately authored review and numerical design check are retained. `inventory.json.gz` contains the complete sequence-only inventory and decompresses to the hash in the manifest. No nonreference measured RNA count is in this directory.

Replay the exact independent sequence-only check from the repository root:

```text
python analysis/programme/retrospective/gill_multiclass_v1/verify_frozen_protocol.py --out work/fresh_gill_multiclass_protocol.json
```

The wrapper verifies all frozen files and reconstructs the original check's local input names in a temporary directory. It requires NumPy but no workbook or network access. It compares its complete result against the retained check. Use a fresh output filename.

Measurement execution remains gated on a frozen selective reader, poisoned-cell fixtures and an implementation commit. RNA1/2 training outputs and their receipt must be committed before RNA3 is decoded. No statistical calibration, biological replicate independence or new gene validation is claimed.
