# B1 independent replay

This folder contains a separately authored numeric reconstruction. It imports no production scientific implementation. The original numerical code was sealed before source/result access; a disclosed schema-only wrapper converts an independently computed exon count into the production summary's exon-ID list for comparison. The failed first attempt is retained.

Actual source replay was authorized by root after the core seal. Do not run another source set, protocol or unselected row through this package. COMT source outcomes remain closed.

Prerequisites: Python 3.12 and openpyxl 3.1.5 (actual version is recorded in the delivery manifest); standard-library exact Fraction arithmetic. No network, Git, author code or workbook macros are executed at replay time.

Use the published Supplementary Data 8 workbook as an external immutable input, the frozen protocol/eligibility and exact primary result. Run from a new output path:

```text
python -B replay_schema_adapter_v1.py --source SOURCE.xlsx --protocol protocol.md --eligibility eligibility.json --result result.json --result-sha256 eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf --authorization actual_authorization.json --output NEW_DIRECTORY
```

Keep this folder's `seal.json` and original files together. The original CLI validates its sealed bytes, protocol, metadata, source/result hashes and scoped authorization before response decoding. Optimized Python is rejected. Existing output directories are refused. The command uses caller paths and does not require the original Windows workspace. The full source workbook and full paper are not duplicated here.

The successful run is `actual_v2`; no additional rerun is necessary for this handoff. The initial `compare_and_run.py` alone intentionally retains its original count/list summary mismatch, so use the disclosed wrapper for complete comparison. Exact source decimals, not openpyxl floating rounding, determine bounds and orders. Positive numeric agreement verifies the defined finite-source functional, not measured biological truth or source error calibration.
