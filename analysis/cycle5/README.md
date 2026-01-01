# Cycle 5: measured joint-function feasibility

This directory records a limited, prospectively specified reanalysis of processed GFP reporter output and relative growth. The scope is an author-preprocessed holdout within one library. It is not independent-study validation or calibrated joint biological inference. No manuscript or submission asset belongs in this directory.

## Authoritative files

`cambray_pilot/protocol_v1.md` is the combined frozen specification. It preserves the reviewed core draft and explicit amendment. Its opening status supersedes historical draft/pending language in retained source files. `cambray_pilot/freeze_v1.json` records the approved inputs, implementation and reviews. Original protocol and implementation bytes are retained for traceability.

The staged implementation is `cambray_pilot_preflight/staged_analysis.py`, with sibling `graph_rules.py` and `masking_v2.py`. It requires CPython 3.11 or later and the standard library; CPython 3.12.14 is the tested version. It rejects optimized Python, changed input hashes and reused output locations. It reads only fields permitted for the requested arm and phase. File hashes and authorization records establish investigator provenance, not a cryptographic security boundary against a deliberate operator.

The source archive is external to Git:

- Source: https://osf.io/download/pm2qa/
- Filename: `Supplementary_Data_15_Integrated_Phenotypic_Measurements.csv.zip`
- Bytes: 130122456
- SHA256: `2758ed63e469c77b42c01e3c022452ef52ebf8d1bd4a86bd4c5b502e44b7bbc0`
- Required member: `Supplementary_Data_15_Integrated_Phenotypic_Measurements.csv`
- Records: 244000

Do not substitute a changed download, different processed mean or missing value. Source data retain their original attribution and terms. The computational license does not relicense third-party source material.

## Reproduction order

Run the synthetic reader suite from a writable copy:

```text
python analysis/cycle5/cambray_pilot_preflight/synthetic_reader_tests.py
python -I -B analysis/cycle5/amendment_diagnostics/test_diagnostics.py
```

The reader tests use temporary invented data and do not open the biological archive. The second suite verifies the added group, target and component summaries. Independent implementation checks and exact arithmetic examples are preserved beside their receipts. Finite checks are not universal proofs.

For biological reproduction, prepare an explicit authorization containing the nine required paths and SHA256 values. Use `portable_handoff/development_predict_authorization.example.json` as a disabled schema example only, replacing its draft protocol hash with the final frozen protocol and its archive placeholder with the authenticated local file. Resolve paths against the checkout root. Execute one phase and one arm at a time with a new output directory:

```text
python analysis/cycle5/cambray_pilot_preflight/staged_analysis.py predict --authorization AUTH_DEV_PREDICT.json --root CHECKOUT_ROOT --arm development --output NEW_PREDICTION_DIRECTORY
python analysis/cycle5/cambray_pilot_preflight/staged_analysis.py evaluate --authorization AUTH_DEV_EVALUATE.json --root CHECKOUT_ROOT --arm development --output NEW_EVALUATION_DIRECTORY
```

The evaluation authorization must pin the prediction seal produced in the first step. That seal pins the numerical predictions and known inputs. Author means, target eligibility and component summaries become visible only in evaluation. If the fixed author-mean check fails, the output is a provenance hold and scientific interpretation stops. Retain that failure.

After a permitted evaluation, the supplemental diagnostic command takes only its exact evaluation bytes and the frozen split:

```text
python -I -B analysis/cycle5/amendment_diagnostics/diagnostics.py --evaluation EVALUATION_JSON --evaluation-sha256 EXACT_EVALUATION_SHA256 --split analysis/cycle5/cambray_pilot_preflight/split_manifest.json --allow-arm development --output NEW_DIAGNOSTIC_DIRECTORY
```

Freeze the complete development report before deciding whether the unchanged question justifies validation. Reserved commands additionally require a pinned development report and explicit unchanged-procedure decision. The full authorization schema and decision fields are in `portable_handoff/README.md`. No reserved missingness or labels may be used to choose whether to validate. An unavailable outcome is not zero function, and a failed prediction screen cannot be retuned on the same reserved data.

## Interpretation and historical replay

The source audit distinguishes nominal preparations from deposited series, one founder from three growth lines, and genotype-aligned traits from paired culture measurements. Twelve cross-series combinations are sensitivity scenarios. Empirical ranges are deterministic sets with uncalibrated biological coverage. Equal group weighting and all 16 group deletions do not create independent experiments. A local square result can be overturned by an observed component bypass and never excludes unmeasured synonymous detours.

The metadata inventory and graph, split, mask, overlap and resolution checks are retained. Some historical preparation scripts resolve fixed sibling paths, use cached author sources or write into their own directory. They must be replayed in a fresh disposable tree with assertions enabled. `portable_handoff/README.md` enumerates these limits. Historical pending labels, absolute source locations and earlier source hashes are records of their original executions, not current commands or claims of a portable raw-assay pipeline. Do not overwrite old receipts to make them appear current.

`limitations_closure_preoutcome.md` and the review folders record the decisions made before new outcomes. Later arm reports and decisions retain failed predictions, unknown states and source holds. The current protocol allows descriptive evidence, an inconclusive test or a provenance failure. Independent biological reproduction and practical small-effect exclusion require additional evidence.
