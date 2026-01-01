# FAS portable replay review

**PASS.** A disposable relocated archive subset with no `.git` directory replayed all intended independent calculations with an empty process `PATH`. The mathematical output equals its original saved JSON exactly, including all 30 exact certificates and 74,124 pair records. The numerical output agrees exactly on all 3,816 selected-field scalar comparisons. Maximum numerical discrepancy is zero. The complete empty-PATH replay took 6.46 seconds on the existing CPython 3.12.14 Windows runtime.

## Changes and retained evidence

The original numerical verifier, mathematical verifier and their original evidence were copied byte-identically into `package/`. No shared source, result or repository file was changed. The mathematical source is entirely unchanged. The numerical offline adapter removes its Git subprocess import, replaces only the Git commit/blob gate with the four fixed scientific/format gate-file hashes, and adds the explicit metadata flag that Git history was not checked. Every scientific helper function is AST-identical. Every statement in the observed numerical comparison after the provenance gate and before the result object is also AST-identical. The only result-object addition is the offline-provenance flag.

`adapter.diff` records the exact change. `preparation_receipt.json` records original/adapted hashes and AST evidence. The wrapper handles paths and gzip decompression, invokes the preserved mathematical verifier and numerical adapter, verifies all 32 archived inputs, and compares their outputs with the original review evidence. It never imports the mutable production analysis modules.

## Adversarial portability checks

Nine checks passed: optimized Python `-O` refusal, `-OO` refusal, overwrite refusal before input access, manifest-tamper refusal, raw-count tamper refusal, adapter-source tamper refusal, compressed-result tamper refusal, missing-input refusal, and a complete relocated replay without Git metadata and with an empty `PATH`. Every failed-input test produced no output directory. The altered disposable inputs were restored to their pinned bytes before the final complete replay. Source inspection found no network-client import or request in either replay path. No claim of operating-system network sandbox enforcement is made.

The wrapper accepts explicit `--repo` and `--out`, works from an unrelated current directory, and uses the current Python executable directly for its subprocesses. Its Windows path normalization includes drive roots and UNC roots. The tested fixture used a local Windows drive, not a live UNC share. Linux/macOS execution was not tested here.

## Handoff

Copy only the contents of `package/` into `analysis/cycle6/result_reviews/`. The relocated test fixture `r/`, decompressed replay outputs and test scripts are review evidence, not required archive payload. The exact command and inventory are in `package/README.md`. Root can retain this report, `adapter.diff`, `preparation_receipt.json`, `portability_receipt.json` and `replayed_empty_path/replay_receipt.json` as supplementary reproducibility evidence.

| Artifact | SHA-256 |
| --- | --- |
| Wrapper `package/replay.py` | `278a6e4424ea35550fed1e137ada1302822b129cf1338dab2307e0ef220603d2` |
| Fixed replay manifest | `2aaf5ceee510b1fb51458e6c733152d1864232ef8ba5062c3791071a1a1982bc` |
| Original numerical verifier | `f5a932166d3cd8ddefc47020e190610b0b845aea9a72a1d6a3531f9693ffd840` |
| Numerical offline adapter | `41cb805ae2cbf0142a35bd7714f347d04e9f71a6eb97a4453cb42e92896c0e39` |
| Original mathematical verifier | `047846795d9383eeef60ad8f74887783f706406fde0f4a75f916bd33af667db2` |
| Adapter diff | `b11dd336e852fabd73f2370748308fa282ef444a5869cc70e19316ef22257df1` |
| Preparation receipt | `7785bdf41f1edd824ba10c4c44af9f3a5eca0cceed57ba0151383686c74dd34f` |
| Portability test receipt | `28d84b34460064f1bd78841ae804c9a0c7c0b5101cf966556422002b7e3575b9` |
| Complete empty-PATH replay receipt | `b5954a3907280c313dfa4f715d5e13909292a447b0b15ee652c9d22b1c604106` |

The original amendment commit remains a provenance identifier. This offline replay deliberately does not assert that it has inspected Git history. No scientific method, endpoint, exclusion, threshold or interpretation was changed.
