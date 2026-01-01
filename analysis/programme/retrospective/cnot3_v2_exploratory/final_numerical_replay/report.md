# Final numerical archive replay review

**Verdict: passed within the unchanged verifier's declared scope.** A clean archive-like directory containing only the 24 pinned computational files reproduced all 59,268 checks and 240,144 numerical comparisons. Its complete verification receipt was byte-identical to the original, with maximum absolute discrepancy `3.397282455352979e-13`. The nine saved models and the negative practical result were unchanged.

The wrapper is standard-library code. It invokes the original independent verifier and reference core without modification, through an explicit argument list and isolated Python mode. The verifier hashes the author implementation but never imports it. There is no Git or network call and no dependency on the original workspace path. All archive file hashes were checked before execution and again after execution.

The source and data manifest includes all six files from the separate cohort review as provenance. Their hashes are verified, while the grouping and eligibility audit is not rerun. This preserves the distinction between numerical replay and source-data validity.

## Observed checks

| Check | Result |
|---|---|
| Clean fixture inventory | Exactly 24 pinned files; no `.git` directory |
| Original verifier and reference core | Copied byte-identically |
| Scientific numerical checks | 59,268 passed |
| Numerical values compared | 240,144 |
| Maximum absolute discrepancy | 3.397282455352979e-13 |
| Prediction / out-of-fold vectors | 9 / 35 |
| Ridge penalty selections | 7 |
| Full evaluation checked | Yes |
| Archived verification receipt | Byte-identical replay |
| Archive input and provenance hashes | All unchanged |
| Standard-library guard tests | 10 tests passed |

The guard tests cover missing inputs, drift in source/data/provenance hashes, changed replay manifest, attempts to reuse an output directory, outputs inside the archive, relative path escapes, duplicate or nonfinite JSON literals, and optimized Python execution. Failed prerequisites do not launch the numerical verifier. The unchanged numerical verifier separately checks its exact runtime and authorization before decoding the scientific input arrays.

## Files and hashes

| File | SHA256 |
|---|---|
| `replay.py` | `71ce3599cc4f9636644eae7bb77841f54fcb0bc987ba5b99e425acf0b747d20b` |
| `manifest.json` | `379b1cdbe4fd394fa6137e45e2fe3eec9b3a0352cda3108d7c0d733ee523a2ff` |
| Unchanged `verify.py` | `896f580c33b8f24e3f44180501c51fd141c2b379b1d6f0631aecf5a75558c246` |
| Unchanged `reference_core.py` | `a5255167a8182d36e107db804d42b87a488e3dc22361cc75c57772a7729c8786` |
| Original and replayed `verification.json` | `9594b34f602488594a855a249a7e34db260de108a7afdf87d1eae3e631840c45` |

The machine-readable handoff receipt records the remaining delivery hashes, fixture inventory and guard-test output. The original verifier's existing prediction tolerances are absolute `1e-9`, relative `1e-10`; its scoring tolerances are absolute `1e-12`, relative `1e-10`. No acceptance tolerances were changed.

## Portability qualification and retained limitations

The verified runtime is Windows, CPython 3.12.14, NumPy 2.5.4, Biopython 1.88. The wrapper itself needs only the standard library, but the existing numerical verifier needs the two pinned scientific packages. Other runtimes are refused, not installed or silently substituted. Other operating systems and BLAS builds remain untested. The wrapper reports the historical maximum difference rather than requiring it to match exactly on every otherwise supported build.

An initial fixture copy encountered the Windows legacy path-length limit. Explicit Windows extended paths resolved it; the final replay ran successfully from a nested archive-like directory. No original computational file changed. An initial manifest placed the cohort-review support files under the results tree; that unexecuted candidate was retained locally and replaced with the root's actual analysis-tree integration before the final replay.

The replay does not reauthenticate historical Git ordering, reproduce the original published cohort, refit the author workflow, recheck experimental units, or establish biological causality. It independently recomputes the existing frozen numerical definitions to verify saved outputs. The negative practical result and unresolved source and measurement limitations remain intact.
