# Independent COMT saved-RF trace audit

Status: **synthetic-ready, actual RF-result access closed pending root authorization**. No actual COMT/source-simulation numeric values or stage outputs have been opened by this reviewer. Public freeze `13b900f704b76b1b3d3059a3bbec8ac1f510b715` is supplied by root, not independently reauthenticated by Git here.

`replay_forest.R` reconstructs source factors and the feature matrix independently using extracted row identities/covariates. It does not import `prep_train_dt`, the author RF adapter or any training/tuning function. It passes the reconstructed matrix to the pinned randomForest implementation and the saved fitted forest. It verifies every retained probability and class, exact record order/identity, source CAF restrictions, every full-row frequency/global-rejection mask, gDNA sample sets, per-sample single-SNP false-call counting and >=4/7 versus >=2/4 rejection thresholds. It checks feature columns, class order, saved seed9 metadata, 500 trees and the saved tuning-to-mtry aggregation.

This is a replay using the same installed prediction engine, not an independently implemented forest engine or independent training/tuning audit. A seed field cannot authenticate the actual RNG history. It cannot establish historical classifier equivalence, biological accuracy or assay calibration. It deliberately performs no new training, threshold selection or subset selection.

The probability tolerance is fixed **before actual access** at absolute 1e-14, permitting decimal serialization of saved probabilities. Class labels and all masks/identities must agree exactly. Feature columns must agree exactly; no dropped/added factor columns are repaired. Full traces are saved before a probability/class failure. Other discrepancies stop without fallback. The original 500-tree classification produces vote fractions; this tolerance does not represent biological uncertainty.

`replay.py` accepts only explicit CLI paths through a trace authorization JSON. It refuses optimized Python, existing output, unapproved actual mode, changed input hashes, broken extraction-to-RF lineage, changed common protocol/code pins, missing/unsealed payloads, changed R executable or package files. It launches fresh `Rscript --vanilla` and retains outputs/log/receipt. Successful root-provided actual authorization must include schema `COMT_RF_trace_authorization_v1`, status `approved_fixed_trace_only`, mode `actual`, `actual_RF_result_access_authorized: true`, `pilot_root`, `verifier_python_sha256` and path/SHA256 inputs for protocol, amendment, integration_addendum, contract, code_manifest, extract_manifest, rf_manifest, runtime_manifest and verifier_r. This is a reproducibility guard, not independent proof of human authorization.

Execution, only after root explicitly supplies and approves actual inputs:

```text
python -B replay.py --authorization AUTHORIZATION.json --output NEW_OUTPUT_DIRECTORY
```

The pinned runtime is R4.6.1/randomForest4.7-1.2 with manifest SHA256 `a97c2057080245e691640ca955af6248aa93102b58341b76f37f473461231c95`. Python3.12 is used for hashes/orchestration; the numerical replay uses base R and randomForest only. There is no installation, network, raw-archive reader, source retrieval or Git operation in the verifier.

The existing wholly invented fixture was copied into the owned `invented_v1` directory. Saved forests were replayed without training. All 440 prediction probabilities and labels matched exactly, all 440 row statuses matched, and both 20-SNP rejection sets/counts matched. This involves 480 invented simulation training rows stored in the earlier models, not actual source simulations. Zero and nonzero probability differences are retained in trace files. Optimization refusal and closed-actual-mode refusal passed. The direct run command was `work/env/Scripts/python.exe -B -X utf8 work/research_cycle8/comt_rf_independent_audit/test_invented.py`.

No failed scientific attempt occurred. Before initial execution, draft syntax was corrected to set print precision through R options rather than a nonexistent write.table argument, and the empty false-call case was explicitly handled. These were unsealed pre-execution edits. No tolerance or numerical method changed after observing actual outputs. The next action is root-supplied actual phase/output paths and explicit access authorization; no actual run is performed by this handoff.
