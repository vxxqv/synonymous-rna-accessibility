# Frozen CNOT3 comparison

The scientific rules are in `protocol.md`; `pre_outcome_freeze.json` pins the reviewed code and computational inputs before the first gene-level outcome decoding. This branch asks whether measured association improves prediction of deposited depletion half-lives relative to the named baseline and composition models. It does not identify a causal pathway rate or provide independent experimental replication.

The data derive from Zhu et al., DOI10.1126/science.adq8587, and the same-source reanalysis by Carneiro and Palhano, DOI10.1002/2211-5463.70277. Source URLs, reference definitions and exact hashes are under `data/programme/cnot3_v1/`. The archived reference FASTA permits exact sequence-group reconstruction without an external sequence lookup. Rebuilding reference features from the genome annotation additionally requires the optional checksum-pinned GTF listed in `provenance.json`; the independently checked extracted feature table is included.

Execution uses CPython3.12.14, NumPy2.5.4 and Biopython1.88. Every scientific phase is launched through `guarded_phase.py`, which checks its source hash and runtime before starting. Outputs must be fresh and outside the archived inputs. Failed runs remain evidence and cannot be silently overwritten or converted into a passing receipt.

The phase order is fixed:

1. `reproduction/reproduce.py`: reconstruct the published Figure2 cohort and contrasts using the archived workbook and `reproduction/authorized_freeze.json`. The historical filename `protocol.draft.json` is the exact immutable operational reproduction definition authorized by that freeze; `protocol.md` is the final scientific protocol.
2. `principal_gate.py`: verify count, arithmetic, ranking, source-field and direction checks. Secondary PTMD warnings remain separate. A failure stops extension for source diagnosis.
3. `cohort/pipeline.py`: apply the fixed original-ID/reference exclusions, exact near-protein groups, resource caps and deterministic split. Require the upstream gate, exact input/code authorization and a complete final receipt. Partial components or a failed coverage receipt cannot authorize training.
4. `prediction/predict.py predict`: fit only on training targets and write all nine model outputs. The strict predictor files contain no depletion target. Root verifies and commits the exact prediction hash before scoring.
5. `prediction/predict.py evaluate`: check the committed prediction and heldout-target hashes before decoding targets. Preserve every comparison, failed criterion, sensitivity and component deletion.

Phase-specific commands and authorizations are retained with completed outputs. Use the exact archived files for replay. An outcome-sensitive modification is a new exploratory analysis, not another attempt at this frozen test. The published reproduction exposes cohort outcomes before modeling, so the result is a locked model-heldout comparison after reproduction, not an untouched experimental validation.

Synthetic checks are executable through `reproduction/test_reproduce.py`, `test_principal_gate.py`, `cohort/test_cohort.py` and `prediction/test_predict.py`. Some fixture runners deliberately require a fresh output location. `verification/` records independently implemented checks, scientific limitations and the printed-unit clarification. The kinetic derivation and exact synthetic checks are under `kinetics/`; their intervals are conditional identified sets, not empirically calibrated confidence intervals.

Only data, computational methods, protocols and verification records belong in this repository. Manuscript sources, publication PDFs and exported figures are delivered separately.
