# Portable replay of the frozen source diagnosis

This wrapper replays the two unchanged independent source-audit scripts using only files already present in the computational repository/archive. It preserves the original5,480-gene cohort,49 ratio-tolerance failures and7,359 finite-IP candidates. It does not revise the failed principal gate or run a model, grouping, new endpoint, alternate formula search or scientific threshold.

Prerequisite: preinstalled **CPython3.12.14**, standard library only. No pip packages, Git, network access, downloads or installation are required at runtime. The wrapper accepts an extracted archive directory without `.git`. It refuses `-O`/`-OO` before arguments or input reads, an existing output path, and output inside the input repository, including resolved aliases.

From any working directory, after root integrates this folder:

```text
python -I -B /path/to/repo/analysis/programme/retrospective/cnot3_v2_exploratory/source_diagnosis/portable_replay/replay.py --repo /path/to/repo --out /path/to/new-external-output
```

The repository/archive must retain the exact paths and bytes in `input_manifest.json`. Eleven files are hash-checked: the two original scripts, source workbook, GeneID-only inventory, original result/exposure/cohort files and four archived evidence files. Seven computational inputs are copied into a fresh internal `stage/decay_cycle4/` sibling layout. The archived expected outputs are **not** copied into the original scripts' output paths. The wrapper first runs `audit_source.py --out .../observed_v2`, then `check_candidate_universe.py`. It retains their stdout/stderr and refuses nonzero exits or nonempty stderr, without retry or fallback.

Complete parsed JSON objects, including all provenance fields, must equal the archived ratio summary, candidate-universe summary and1,879-record exclusion ledger. The mismatch CSV must match every header, field and all49 data rows. No fields are omitted from comparison. JSON is compared semantically because platform line endings need not be byte-identical; byte identity is reported separately. All repository and staged immutable-input hashes are rechecked after execution. A failure emits `FAILURE.json` and no success receipt.

On success, the fresh external output contains `receipt.json`, both original process logs and their reconstructible staged inputs/results. The receipt states the original script hashes, exact comparisons and preserved failure. The two scripts' source hashes remain `fa3dfe63521b26fc23a3c8850ab2b5e2d5f99be47a23a6d030a3ce3858256be0` and `21f02e6e2fce5e4a6ca39e2ff1dafb58a3a20e2f70de1d449f4bdab2af0f576a`.

One isolated archive-like replay passed locally on Windows with only the eleven pinned input files and no Git metadata. Both complete summaries and the exclusion ledger also happened to be byte-identical on that host; all CSV content matched. Overwrite, in-repository output and both optimization flags were refused. `verification_receipt.json` records these checks, and `verified_replay_receipt.json` is the complete successful wrapper receipt. Linux execution is not claimed by this Windows receipt. No original audit source, output, repository input or scientific protocol was edited.

The result remains source diagnosis: all5,480 ratios are numerically compatible with the previously examined10-significant-digit intermediate, but the generating procedure is unauthenticated and the published5,481 count remains unresolved. The replay verifies this archived evidence; it neither closes that source gate nor authenticates biological measurements.

Root should integrate **only** the six files listed in `integration_manifest.json`. Local `archive_fixture/`, `verification_archive/`, `verification_run/`, fixture-builder and manifest-builder files are development evidence, not release inputs. This prevents duplicate source workbooks or temporary staged outputs entering the archive.
