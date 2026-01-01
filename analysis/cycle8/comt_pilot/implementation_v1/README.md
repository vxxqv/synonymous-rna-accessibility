# COMT A1 staged implementation, pending final independent review

No real count, simulation-summary value or phenotype outcome has been decoded during implementation. Metadata-only contract construction and invented fixtures are permitted; real execution still requires the parent to approve and freeze the exact protocol, adapter and authorization files.

The implementation calls the held model kernel rather than copying or revising it. It also calls unchanged selected author RF helpers in the installed task-local runtime. Historical classifier-state equivalence is not claimed. See integration_addendum.md for exact scope, denominator and coupling limitations.

Each stage accepts:

```text
python driver.py --phase extract --authorization extract_authorization.json --out NEW_DIRECTORY
python driver.py --phase rf --authorization rf_authorization.json --out NEW_DIRECTORY
python driver.py --phase normalize --authorization normalize_authorization.json --out NEW_DIRECTORY
python driver.py --phase predict --authorization predict_authorization.json --out NEW_DIRECTORY
python driver.py --phase evaluate --authorization evaluate_authorization.json --out NEW_DIRECTORY
```

Authorization is a JSON object with schema COMT_A1_authorization_v1, the exact phase, status approved_frozen_preoutcome, historical_reproduction_claim false, outcome_access_authorized true, and inputs mapping each role to an explicit path and SHA256. These fields are a reproducibility gate, not a substitute for the parent's human-authorized freeze decision. Every phase pins protocol, amendment, integration_addendum, contract and code_manifest. Extraction also pins archive and runtime_manifest. RF pins extract_manifest and runtime_manifest. Normalization pins extract_manifest and rf_manifest. Prediction pins normalization_manifest. Evaluation pins normalization_manifest and prediction_manifest and requires sealed_prediction_sha256 equal to the authorized prediction manifest hash.

The code manifest uses schema COMT_A1_code_freeze_v1, files mapping relative paths from this implementation directory to SHA256, and named_input_hashes for protocol, amendment, integration_addendum and contract. It must include every Python module, apply_source_rf.R, ../proxy_models.py, ../source_helpers.R and both held author helper files. The caller creates this exact freeze only after review; there is deliberately no self-authorizing real-data example.

Every upstream stage must match the same five common protocol, amendment, integration, contract and code pins. The RF stage additionally matches the extraction runtime pin. Normalization requires the RF stage to reference the identical extraction seal. All consumed files must be listed in their prior seals, so an unsealed payload cannot be used. Parsing failures retain the current member and completed-member exposure ledger, without claiming that an interrupted member's numeric exposure was zero.

Input-role paths are explicit. The implementation resides immediately below the reviewed pilot directory so its held kernel, helpers and task-local runtime can be located without a network lookup. Runtime reconstruction from a clean machine requires the pinned R binary packages or another separately reviewed exact environment; this is not an all-stdlib, offline scientific replay claim. The source archive is not redistributed by this package. R is invoked with --vanilla and no package installation occurs during scientific stages.

Run synthetic tests with unittest discovery on test_implementation.py. verify_RF_adapter.py creates a new invented simulation/read fixture and exercises the RF adapter, preserving its receipt in synthetic_RF_v1; it refuses to overwrite that directory. The inherited proxy kernel has its separately held 550-check receipt. These checks test implementation behavior and do not constitute evidence for biological assay substitution.

verify_staged_cli.py writes a disposable synthetic_stages_v2 directory and exercises sealed prediction then evaluation using 125 invented variants at 53 sites. Its earlier v1 receipt is retained as history from before the strengthened provenance guards. verify_source_normalization.R writes an invented two-allele source-helper result; verify_source_normalization.py compares this to independent exact rational arithmetic. These scripts never accept the biological archive as a test fixture.

All 125 targets, failed fits, missing preparations, denominators, raw representation mappings, record filters, model coefficients, ranks, condition diagnostics, out-of-training-range flags, source preparation vectors and same-cohort prediction errors must remain available. The reader rejects unexpected schemas and uses no fallback to published scores. Failed source parsing or fitting must not trigger ad hoc subsets, another seed or another source filter.
