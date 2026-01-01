# Cycle9 offline numerical replay

This bounded computational package replays three previously completed independent numerical checks. It does not execute the author prediction pipeline, source classifier, new biological outcomes or acquisition. CPython3.12 and its standard library are required; the tested patch version is3.12.14. No Git, network, R, numerical optimizer or additional package is needed. The current biological interpretation and all earlier failed criteria remain unchanged.

The archive root contains `analysis/cycle9/FILES.json`, this interface, and a strictly enumerated payload. Copying that subtree into a clean computational checkout preserves every archived byte. The caller must supply the trusted manifest SHA256 from the reviewed staging receipt. The exact command, with an output directory that does not exist, is:

```text
python -B analysis/cycle9/replay_cycle9.py --repo ARCHIVE_ROOT --out NEW_REPLAY_DIRECTORY --manifest-sha256 TRUSTED_FILES_JSON_SHA256
```

For a byte-only check add `--verify-only`. The output must be outside the archive and its ancestors. `-O`, `-OO`, an existing output directory, modified source bytes, an unrecognized path, a missing input or an uninventoried Cycle9 file is rejected. Failure after materialization leaves a failure receipt and child logs. Never reuse that output as a fresh run.

The interface verifies stored hashes and sizes, decompresses the explicitly marked gzip entries, verifies logical byte hashes, and restores a minimal logical workspace. It then executes four unchanged independent scripts: the cPRIME persistence check, fresh COMT local-order check, control source-cell reconstruction and its additional common-C joint-link check. All parsed scientific receipt fields must match the archived reference objects, not just a success flag. Absolute input path labels are canonicalized to their `work/` suffix; no hash value or scientific field is dropped. All60 cPRIME profile CSV rows must also agree.

Expected coverage is396 paired rows,60 profiles,720 assignments and2160 folds for cPRIME;134 pairs and17152 joint cases for COMT;18 numeric source cells,27 ratios,162 quotients and729 shared-C completions for the three-state control. The control's extra checker compares1458 joint quotient links. These are computations on dependent observations, not independent biological replications.

Only four Cycle8 files are included as dependencies: B1's complete selected result, B2's assignment result, the existing COMT source-normalized assay summaries and the fixed COMT source/target contract. Their original and logical hashes remain recorded. Running the full Cycle8 materializer would add irrelevant files and is unnecessary here. The only added primary-data workbook is the explicitly whitelisted Supplementary Data6; its unchanged bytes and source attribution are recorded. The replay reader permits only the already authorized18 numeric cells. No other source workbook, paper, figure, renderer, manuscript or publishing asset is in this package.

## Historical files and the cPRIME guard amendment

Original protocols, source, results, review receipts, guard failures and manifest-correction history are preserved. Historical reference receipts are archived but are not restored over the fresh checker output paths. The manifest records every restoration mapping.

The original cPRIME independent checker pins the pre-guard production source `f99e6d1c...`. Current production source `379e6647...` adds an optimization refusal. The exact historical bytes were recovered by deleting only that documented guard and verifying the original full SHA256. Both versions are archived. The historical version is restored only as evidence for the unchanged independent checker's hash comparison; neither production source is executed. Current guarded production bytes are separately restored under `guarded_current/`. This is an explicit historical compatibility mapping, not an altered freeze or a claim that the current production source remained unchanged. The replay interface and every executed checker reject optimized Python.

Archived authorizations and receipts may contain their original local path labels and earlier pending statuses. They document history; the replay does not rewrite them or claim to reproduce historical timing, access order, user authorization or the originating authors' environment. Current replay status is the new `replay_receipt.json`. The old control execution report's pending independent-review line is superseded by the archived root source-cell and joint-link PASS receipts. No classifier or source-normalization reconstruction is performed.

Tested here on Windows with exact bytes; Linux/macOS require the same standard-library interface and must produce equal parsed scientific objects. Platform-specific line endings and absolute path labels are not falsely claimed byte-identical. This package contains the needed inputs for these finite numerical verifications only; it is not an all-offline reproduction of the full originating studies or the entire research programme.
