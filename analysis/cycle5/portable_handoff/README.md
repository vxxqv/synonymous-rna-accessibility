# Cambray bounded computational handoff

This handoff is documentation and synthetic verification. It authorizes no biological execution. The reviewed candidate modules are unchanged. The protocol entry is a snapshot of the current draft, not a final freeze. Root must finalize the protocol and authorization after independent review. The downloadable source contains measurements, so do not open it outside the staged authorization.

## External source and executable requirements

The author download is https://osf.io/download/pm2qa/ . The deposited filename is `Supplementary_Data_15_Integrated_Phenotypic_Measurements.csv.zip`, containing `Supplementary_Data_15_Integrated_Phenotypic_Measurements.csv`. Its authenticated size is 130,122,456 bytes and SHA256 is `2758ed63e469c77b42c01e3c022452ef52ebf8d1bd4a86bd4c5b502e44b7bbc0`. It contains 244,000 records. These facts are copied from the archived extraction and provenance receipts; the handoff does not open the ZIP. Use the stable OSF link, not the expired signed redirect recorded during download. A nonmatching download must stop execution.

The staged runner and its 13-test suite use only the Python standard library. They require Python 3.11 or later because they use `hashlib.file_digest`; the tested interpreter is CPython 3.12.14, 64-bit Windows. The exact interpreter/platform and copied-source hashes are in `synthetic_candidate_receipt.json`. Cross-platform numerical byte identity has not been demonstrated. Prefer the tested version for final reproduction, use normal optimization, and never set `PYTHONOPTIMIZE` or use `-O` or `-OO`. The runner refuses these modes before imports, file reads or output creation.

Historical full metadata checks additionally use Biopython 1.88; the alternate split checker and original synthetic linear-algebra test use NumPy 2.5.4. Those are the installed versions observed for this handoff, not new installations. The staged pilot itself does not require either package. No Git or network connection is used by the staged runner once all exact local inputs are available. This does not make the whole project an offline reconstruction of the author's raw assay processing.

## Portable pilot inventory and paths

Keep `staged_analysis.py`, `graph_rules.py`, `masking_v2.py` and `synthetic_reader_tests.py` together. Suggested repository location: `analysis/cycle5/cambray_pilot_preflight/`. Preserve the frozen `split_manifest.json` and `prediction_mask_v2.json` there. Keep `components.jsonl` and `synonymous_edges.tsv` under `analysis/cycle5/cambray_eligibility/`, and the finalized protocol under `analysis/cycle5/cambray_pilot/`. The metadata and code hashes are in `manifest.json`.

Only nine input roles are required to execute the pilot. Every role is an object with exactly the operational fields `path` and `sha256`; paths are absolute or resolved relative to the explicit `--root` argument. The three code roles must resolve to the actual executing files. The archive can remain outside the Git checkout and be named by an absolute local path.

| Role | Input |
|---|---|
| `protocol` | Root-frozen protocol bytes |
| `split` | `split_manifest.json` |
| `mask` | `prediction_mask_v2.json` |
| `components` | Sequence-only `components.jsonl` |
| `edges` | Sequence-only `synonymous_edges.tsv` |
| `archive` | Exact external author ZIP |
| `graph_rules` | Executing `graph_rules.py` |
| `masking` | Executing `masking_v2.py` |
| `runner` | Executing `staged_analysis.py` |

`development_predict_authorization.example.json` is an exact schema example with real candidate hashes, a placeholder archive path and `execution_authorized: false`. It cannot authorize an outcome read. Root must set the final protocol hash, actual paths and approved phase/arm after review. Top-level fields are `schema: "Cambray_staged_v1"`, `execution_authorized: true`, `allowed_phases`, `allowed_arms`, `expected_records: 244000`, and `inputs`. Use narrow allowed lists for each invocation.

The synthetic suite was replayed from fresh byte-identical module copies under `isolated_synthetic_replay/`, with invented data only. Run an integrated copy as:

```text
python analysis/cycle5/cambray_pilot_preflight/synthetic_reader_tests.py
```

The suite creates and removes synthetic temporary directories next to its source; use a writable temporary copy when the archive is read-only. It never opens the biological ZIP.

## Authorized execution sequence

All commands below are examples for a future explicit root authorization, not commands executed in this handoff. `REPOSITORY_ROOT` is the extracted computation repository. Output directories must not already exist. Relative authorization paths resolve relative to the current directory; paths inside authorization resolve relative to `--root`.

1. Development known-corner extraction and predictions:

```text
python analysis/cycle5/cambray_pilot_preflight/staged_analysis.py predict --authorization AUTH_DEV_PREDICT.json --root REPOSITORY_ROOT --arm development --output NEW_DEV_PREDICTION_DIRECTORY
```

This reads only the seven allowed replicate fields at globally unmasked known corners. It writes `known_inputs.json`, `predictions.json` and `seal.json`. Preserve all three. Root records the seal hash before proceeding. The outcome target's finiteness, author means and component summaries are not read in this phase.

2. Create a new authorization permitting development `evaluate`, retaining the same nine hashes and adding `inputs.prediction_seal` with the explicit path and SHA256 of that development `seal.json`.

```text
python analysis/cycle5/cambray_pilot_preflight/staged_analysis.py evaluate --authorization AUTH_DEV_EVALUATE.json --root REPOSITORY_ROOT --arm development --output NEW_DEV_EVALUATION_DIRECTORY
```

This verifies the seal and both sealed numerical files before reading the required arm outcomes. Missing target data are retained as failures. A larger-than-frozen author-mean discrepancy or missing mean with complete preparations creates `held_for_author_mean_provenance`; no interpretation or validation decision should be generated from a held arm. A completed arm writes `unsealed_inputs.json`, `evaluation.json` and `evaluation_receipt.json`. A held arm writes the first two only.

3. Root freezes the development report and decides whether an unchanged validation procedure is justified. A reserved-arm authorization additionally pins `development_report` and `development_decision`. The decision must contain:

```json
{
  "status": "development_report_frozen",
  "unchanged_validation_procedure": true,
  "protocol_sha256": "EXACT_FINAL_PROTOCOL_SHA256",
  "report_sha256": "EXACT_DEVELOPMENT_REPORT_SHA256",
  "code_hashes": {
    "graph_rules": "EXACT_GRAPH_RULES_SHA256",
    "masking": "EXACT_MASKING_SHA256",
    "runner": "EXACT_RUNNER_SHA256"
  }
}
```

This is an explicit investigator decision. The runner verifies its hashes and status, but cannot determine from prose whether the scientific decision was justified. It is not a cryptographic approval identity system.

4. Use the same predictor on the reserved arm, without new models or thresholds:

```text
python analysis/cycle5/cambray_pilot_preflight/staged_analysis.py predict --authorization AUTH_VALIDATION_PREDICT.json --root REPOSITORY_ROOT --arm reserved_validation --output NEW_VALIDATION_PREDICTION_DIRECTORY
```

5. Preserve the reserved prediction seal and add its path/hash to a new evaluation authorization, keeping the development decision/report unchanged:

```text
python analysis/cycle5/cambray_pilot_preflight/staged_analysis.py evaluate --authorization AUTH_VALIDATION_EVALUATE.json --root REPOSITORY_ROOT --arm reserved_validation --output NEW_VALIDATION_EVALUATION_DIRECTORY
```

Archive all authorizations, both stage outputs in each arm, final protocol/source hashes, explicit development decision and independent review receipts. Hash verification reads opaque archive bytes; it does not decode phenotype values. Evaluation output remains conditional on processed reporter/growth measurements and the declared observation graph. It is not raw replicate replay, full-genotype-space accessibility or independent biological replication.

## Metadata replay and historical limits

The pilot can use the reviewed hash-pinned metadata without running historical preparation scripts. To independently reconstruct metadata from the archived sequence inventory, retain the 30,751,838-byte `sequence_metadata.tsv`, `annotated_squares.json`, `binary_faces.jsonl`, `components.jsonl`, `synonymous_edges.tsv`, `previously_exposed_ids.json`, all split/mask manifests, historical code and synthetic receipts. Their exact inventory and hashes are in `manifest.json`. Sequence inventories and exposure IDs contain no measured phenotype values.

Preserve sibling folder names `cambray_eligibility`, `cambray_pilot_preflight` and `cambray_pilot` in a disposable replay tree. Several historical scripts write into their own folders, overwrite receipts or refuse existing outputs. Do not run them in the immutable archive. Their use of `assert` is not protected by a module-level optimization refusal; execute them only through a guarded replay driver or an interpreter check that aborts when `sys.flags.optimize != 0`. The authoritative staged runner already has its own refusal. Do not silently alter historical source hashes to retrofit guards.

Specific portability gaps:

- `metadata_split.py` resolves `../cambray_eligibility/` and `../cambray_pilot/protocol_draft.md`. It stores the protocol draft hash. A later protocol freeze will change that recorded historical hash without changing seed or corner assignments. Preserve the historical manifest and compare scientific unit content explicitly instead of pretending a new full-file hash reproduces it.
- `verify_split.py` also expects the sibling eligibility folder and writes `split_corroboration.json`; it needs NumPy and Biopython. `test_synthetic.py` writes its original receipt; `test_masking_v2.py` and `masking_v2.py` require historical manifest files and fresh destinations. These are reconstruction scripts, not arbitrary-path production CLIs.
- `audit_graph.py` additionally expects `HERE.parents[1]/literature_expansion/cambray_generate_seed.py` and reads the author's literal standard-code table. Its authenticated source is https://osf.io/download/5y7g9/ with SHA256 `17f201c7048e4570e8d0071f9ab8047230e5e7027f8be43eabb9426f41e393d2`. That is not present in the minimal pilot inventory. `check_graph.py` and `square_design_check.py` can corroborate existing inventory from their sibling files, but do not replace all author-context provenance checks.
- `recover_sequences.py` expects the primary archive under its own `sources/` folder and a prior `literature_expansion/cambray_zip_head.bin`. It can download the primary archive if absent and reconstruct the historical 503-row exposure boundary. That compressed prefix contains old exposed measurements and is deliberately not required or copied by this handoff. Preserve its hash and `previously_exposed_ids.json`; do not claim from-scratch exposure reconstruction without obtaining the prefix and authorization. The current pilot conservatively excludes all seed1 designs regardless.
- `finalize_receipt.py` depends on cached author code, OSF inventories and several source files under fixed `literature_expansion` and `sources` paths. It is a historical provenance assembler, not a portable raw-data pipeline. Its previously saved receipts and source URLs are retained.
- The historical v1 README includes older pending decisions and the phrase “separate-person review”. It is preserved byte-for-byte as history. Current review means separately authored analyses with shared project context, not blinded external human review. The final protocol and this handoff supersede those status statements.

Do not redistribute full primary articles or a publisher challenge page as article text. Retain primary DOIs, stable author-deposit URLs, content hashes where authenticated and explicit unavailable-source statements. Keep the 130 MB phenotype ZIP external to Git. Public URLs may later change or become unavailable, so source retrieval is a remaining external dependency unless a matching local copy is supplied. No new download or biological reading occurred for this handoff.
