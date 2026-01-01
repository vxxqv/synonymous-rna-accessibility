# Offline independent FAS replay

Copy this directory's contents to `analysis/cycle6/result_reviews/` in the computational archive. Preserve all files byte for byte. The replay requires only Python's standard library, the 12 archived source count files, the compressed result, the frozen scientific/format inputs, and the original independent review evidence listed in `replay_manifest.json`.

From any working directory, run:

```text
python /path/to/checkout/analysis/cycle6/result_reviews/replay.py --repo /path/to/checkout --out /path/to/new-audit-directory
```

The output directory must not already exist and must be outside the computational archive. No Git installation, `.git` directory, network connection, scientific Python package, author source module import, refitting or optimization is required. Python 3.10 or newer is accepted; the recorded complete replay used CPython 3.12.14 on Windows. Other operating systems were not exercised in this portability review. Do not use Python's `-O` or `-OO` options. Windows drive and UNC roots are normalized to extended paths by the wrapper.

The wrapper verifies the fixed manifest, all 32 input byte hashes, the compressed result, and the exact 29,145,849-byte decompressed result before running either verifier. It preserves the original independent verifiers and uses a numerical adapter whose only functional change replaces the Git-history gate with fixed archived-byte hashes. An explicit `git_history_checked_in_this_replay: false` field prevents conflating offline numerical replay with authentication of the historical commit. The original scientific and format freezes remain unchanged.

The replay runs all ten numerical synthetic safeguards, the exact raw-count/normalization comparison, and the complete mathematical verification of 30 certificates and 74,124 endpoint-pair records. It then compares the reconstructed evidence with the original saved independent evidence. Mathematical JSON must match completely. Numerical JSON must match after removing only runtime, verifier-source hash and the additional explicit offline-provenance flag. Input files are hash-checked again at the end.

The output retains the decompressed result, `synthetic.json`, `numerical.json`, `mathematical/verification.json` and `replay_receipt.json`. Timings and their enclosing receipt hashes vary by run; mathematical content and the declared exact comparisons must not. A partial output from a failed invocation is retained for inspection and cannot be overwritten by a later invocation.

## Minimal archive inventory

| Path | Required contents |
| --- | --- |
| `analysis/cycle6/result_reviews/` | This wrapper, fixed manifest and both reviewer subdirectories |
| `analysis/cycle6/result_reviews/numerical/` | Exact original `verify_counts.py`, derived `verify_counts_offline.py`, original observed and synthetic receipts |
| `analysis/cycle6/result_reviews/mathematical/` | Exact original `verify_result.py` and original `verification.json` |
| `analysis/cycle6/fas_pilot/` | Nine scientific-freeze inputs plus `freeze.json`, `format_freeze.json`, `run_headerless.py` and `format_amendment_v1.md` |
| `data/cycle6/fas_counts/` | All 12 exact original count-file blobs |
| `results/cycle6/full_result_v1.json.gz` | Fixed compressed result |

The manifest gives every exact path, size and SHA-256. This README and the wrapper are additional inventory items; the wrapper pins the manifest hash. No manuscript, figure, publication package or downloaded article is required.

## Interpretation limits

Successful replay authenticates these fixed computational inputs and their arithmetic. It does not establish independent biological validation, calibrated sampling uncertainty, source laboratory provenance beyond the archived metadata, or causal interactions. The 27 ratio scenarios share nine outputs and three technical inputs. Their envelope is not a confidence interval. The secondary point is a corrected headerless source-formula reconstruction, not an authenticated historical author result.

The Git-dependent original numerical verifier is deliberately retained as historical evidence. Use `replay.py` for the offline route. The original Git-history check remains available only when its corresponding commit exists locally.
