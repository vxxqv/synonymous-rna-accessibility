# Final Cycle8 archive guard review

Verdict: PASS for the reviewed final inventory, path validation and synthetic guard scope. No remaining portability/dependency blocker was found. No COMT or other scientific replay was run for this review. The earlier unsafe-copytree snapshot and its findings remain preserved in the parent review directory.

## Reviewed immutable points

| File | SHA256 |
|---|---|
| materialize_cycle8.py |37b2199996d90735cfe15b109ee24c483358e6ce5f8977d5afc4dfd08b6f9041|
| verify_materialized_cycle8.py |eb77d6732eeca94900571545086eb40e2765dab4b3972d64daa85b4fcc0b964a|
| REPLAY_FILE_INVENTORY.json |1826554f35e04b9d5b4761992f0aa94e3c23ab7129438b079fb8557018b9a8df|
| README.md |699674c1f22f072d7d59df1fdae0c1e2f6b9855934ac347d9eba01998a9759dc|
| CURRENT_VERIFIED_STATUS.md |aa2b29d38cbca85c3c87f1cd2e253e808ffe2b5d34637c6ad8767903cf12cda4|

The authoritative inventory covers725 files and49,997,497 stored bytes. Every file's hash and byte length matches. The only physical non-bytecode file omitted from its own list is REPLAY_FILE_INVENTORY.json itself, avoiding a self-referential hash. No symlink or normalized duplicate was encountered in the listed files. All listed paths are relative and confined beneath analysis/cycle8.

The lossless logical manifests now contain282 entries:90 A,99 BDE and93 CFG. Their two-entry reduction is precisely the removal of the two unnecessary whole notebooks from unpublished archive copies. All282 stored and logical hashes/lengths pass. The five executed scientific checker paths are among those logical entries, so the combined verifier rechecks their hashes before execution. The complete725-file inventory also seals the earlier computational source/protocol/history files outside those282 logical-result entries.

No .ipynb,.pyc,PDF,image,TeX,DOCX or PPTX is listed in the inventory. Execution logs, correction diffs, classifier RDS objects, selected numeric records and the original cPRIME workbook are intentional computational inputs/evidence. The materializer copies only listed files; ignored local caches and other stray source-tree files cannot enter by recursive directory copying. The original source notebooks and their hashes/extracted code remain provenance; no notebook cell or output was opened in this audit.

## Independently tested guards

check_guards.py performed34 invented checks without scientific inputs. Direct safe_join tests reject empty/nonstring,absolute,UNC,drive-qualified,backslash,colon,current-directory,parent traversal,repeated-separator and trailing-separator paths. A nested ordinary path is accepted.

Tiny invented repositories exercise the actual materializer CLI and reject duplicate inventory paths, Windows case aliases, absolute/parent inventory paths, altered source hashes, duplicate logical destinations/case aliases, logical/archive traversal, unknown compression and logical files missing from the sealed inventory. All these negative fixtures fail before creating their requested output directory. Existing output and output-inside-repository guards also pass. Both materializer and combined verifier refuse -O and -OO before any argument/input access. The ancestor-output fixture is already existing and therefore is rejected by the earlier no-overwrite guard; the explicit ancestor-containment test remains additionally evident in the code.

Actual inside/outside directory symlinks were attempted only within the owned synthetic fixture area. Windows returned privilege error1314 for both creations. Their rejection is therefore not claimed as an executed real-symlink check. Static review confirms safe_join resolves containment and inspects every path ancestor for symlinks; the strict materializer inventory is also confined. No cross-platform filesystem test or concurrency/TOCTOU proof is claimed.

The main filesystem test script was invoked once. Importing the reviewed materializer created a disposable ignored Python cache; that test-created cache was removed afterward, with no source/result byte change. Future invocations should use python -B. The synthetic fixtures and receipt are preserved under tests_v1. No original scientific suite, classifier training, normalization or prediction was executed.

## Dependencies and current documentation

requirements-tested.txt now includes scipy==1.18.1, alongside numpy==2.5.4,sympy==1.14.0 and openpyxl==3.1.5. The README correctly requires CPython3.12.14 and distinguishes exact certificate checking from the frozen SciPy solving implementations and separately pinned R classifier runtime. It documents a fresh materialized workspace plus a fresh verification directory, preserving D2's sibling-output convention. It explicitly limits COMT replay to frozen extracted counts/masks and denies recovery of historical classifier training/raw source authentication.

The COMT criterion is correctly reported as failed and numerical reconciliation as passed. Current-result indexing supersedes historical pre-outcome language without changing frozen protocol bytes. The documentation retains descriptive uncertainty, selected-source ascertainment, the marker-domain correction, row/normalization dependence, finite FAS support scope and conditional51/805 kinetic classification. Nothing claims external biological replication or converts missing source files into a biological null.

Dao source_license.md now identifies the exact journal/DOI,author commit,MIT statement and preserved source attribution. The local DOI-matched primary XML's Data availability and Code availability text explicitly states MIT for the linked data and analysis scripts/notebooks. The notice honestly records that no separate LICENSE file was found and invents no copyright holder. Root THIRD_PARTY_NOTICES.md now distinguishes all three Cycle8 sources from the repository's new-code MIT notice. This closes the source-record gap identified in the earlier review; no legal certification beyond checking those source statements is asserted.

## Remaining interpretation and assurance limits

This is verification of byte-preserved packaging and exercised guards, not a new scientific replication. The complete numerical/symbolic replays were already performed by root and reviewed via their saved receipts; they were not repeated here. The archive itself still needs its final top-level checksum/version/commit record. The materializer's in-archive inventory is a content-consistency mechanism, not a cryptographic signature against an adversary rewriting both code and manifests. Standard source/archive authentication remains separate. No remaining issue found in this bounded review requires changing scientific functions, thresholds, predictions or failed-result history.
