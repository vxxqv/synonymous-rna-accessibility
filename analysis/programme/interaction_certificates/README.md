# Portable A3 interaction-certificate evidence

This package replays the frozen author computations, the separately authored finite reference, and the final revision-2 accounting diagnostics. It needs Python 3.12 with its standard library and no Git, network, scientific packages, original workspace or datasets. The successful packaged replay used CPython 3.12.14. Other Python versions have not been certified for the source-line accounting diagnostics.

Copy these four files together into any directory: `replay.py`, `computational_sources.zip`, `integration_manifest.json` and this `README.md`. From there run:

```text
python -B replay.py --out NEW_DIRECTORY
```

`NEW_DIRECTORY` must not exist. The entrypoint refuses overwrite and optimized Python execution (`-O` or `-OO`), because the frozen computations use assertions. Expected runtime is seconds to about one minute on a normal desktop. There is no Internet or Git-history query. A failure is retained in the new output directory and returns a nonzero exit code. A passing run writes `receipt.json`, with exact comparisons and input hashes before and after. It does not modify the four package files or any original project files.

## What is preserved and replayed

`computational_sources.zip` contains the original source bytes, contracts, proofs, prior-work notes, freezes, finite receipts, independent review and revision history. Every member is listed and hashed in `integration_manifest.json`. The imported additive dependency is pinned to SHA-256 `96c6d35a6243dcf1ba0dbd08100d82dbbf553a583e3844e9ac6a31646d5ec63a`. No primary literature PDF, manuscript, figure, publication export, expression data, Git directory or bytecode cache is included. The review's small PDF retrieval receipt is metadata only; the PDF itself is excluded.

The new output has two separate trees:

- `preserved/` retains every archived file byte-for-byte, including all historical results.
- `reconstruction/` receives only static inputs first, in the relative layout expected by the unchanged source. Fresh author, reference and revision-2 outputs are then created there. Historical results are never copied into active output paths before a run.

The adapter calls the unchanged creator functions in their original order. It separately executes `independent_review/reference.py`, whose mathematical oracle uses exhaustive max-min Floyd closure rather than the creator solver. The reference imports creator code only for its explicitly labelled storage/contract diagnostics. Finally the adapter imports and calls `revision_2/targeted_regression_v2.py:suite`, leaving that source's historical Git-dependent `main()` unexecuted. This is a filesystem and provenance adaptation, not a change of mathematical methods. Historical freeze keys with Windows separators are interpreted as relative paths on any host without rewriting the frozen JSON files.

All scientific fields are compared exactly, including counts, statuses, qualifications and list order. Only the creator receipt's top-level `runtime_seconds` and `software`, and the corrected diagnostic receipt's top-level `runtime_seconds` and `provenance`, are excluded. The independent reference output is compared in full. The creator's original full replay core is also compared exactly. Historical receipts retain their original software, commands, paths and Git metadata as records; those historical paths are never followed. A portable run explicitly reports `historical_git_verified: false`. Hash equality preserves that evidence but does not newly authenticate historical repository events.

## Current proof and correction history

The current mathematical and implementation contract is in `preserved/work/interaction_cycle3/revision_2/claim_and_proof.md` and `complexity_contract.md`. `correction_ledger.md` and `trace_erratum.md` are part of this record. The original proof and review are historical evidence, not the current total-cost claim. Earlier review statuses that require a storage qualification remain unchanged rather than being silently rewritten.

Revision 2 adds full-coordinate witness storage, all retained envelope outputs, natural-order width, coordinate bookkeeping, caller component restriction, exact rational bit-cost qualifications and explicit exhaustive-branch costs. In particular, the original implementation does not automatically trim baseline components and does not accept an elimination-order argument. The large implicit product is a representation example, not a general scalability benchmark.

The initial diagnostic `revision_2/results.json` has inflated active-list counters from repeated Python line events. It is preserved only as history. `targeted_regression_v2.py`, `freeze_v2.json`, `results_v2.json` and `manifest_v2.json` are the corrected accounting evidence. The correction changes neither the mathematical outcomes nor the held creator source.

The author suite has 2,818 parameter inputs and 75 small MaxCut source graphs. The reference has 1,296 Q2 parameter inputs but 929 distinct fields, and 32 deterministic Q3 parameter inputs but only seven distinct fields. It records 76,544 projection checks, nine additional one-factor checks and 217 scaled switch checks. These are finite corroboration of the proofs. Generic projection certificates and max-sum variable elimination are established methods. Sparse-interaction hardness remains conditional on the cited cubic-MaxCut theorem, whose full proof was not independently re-proved here. No new general algorithm, priority, calibrated uncertainty or biological result is claimed.

## Integration

Archive the four distribution files together. A retained successful `receipt.json` can accompany them as execution evidence, but its generated `preserved/` and `reconstruction/` duplicates need not be archived again. The integration manifest identifies every source and historical receipt and states which corrected documents are authoritative. No mathematical source editing is needed to integrate or replay this package.
