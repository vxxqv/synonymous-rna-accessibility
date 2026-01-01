# Portable robustness certificate replay

The original `verify_result.py` and its receipt remain unchanged. Use the derived `verify_portable.py` for archive and CI replay:

```text
python analysis/cycle7/robustness_result_review/verify_portable.py --repo /path/to/checkout --out /path/to/new-review-directory
```

`--root` is an alias for `--repo`. The output directory must be new. By default the verifier reads `results/cycle6/full_result_v1.json.gz` and `results/cycle7/robustness.json.gz` relative to that checkout. An explicit `--result /path/to/robustness.json` or `.json.gz` overrides the successor path. The exact decoded successor bytes remain pinned to SHA-256 `861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa`. The compressed container's actual hash is recorded and checked again after the run; different valid gzip metadata does not change the authenticated decoded result.

The archive must also contain `analysis/cycle7/freeze.json`, all six source inputs it lists, and the Cycle 6 ordinal core and sequence contract. The verifier reads these sources only to validate their hashes. It does not import them. No Git, network, NumPy, SciPy or other scientific package is needed. The implementation uses standard-library exact rational arithmetic. It rechecks the existing certificates and does not rerun the optimizer or calculate a new sensitivity analysis.

Only `main` differs from the original verifier. The changes add the root alias, default successor path, gzip decoding and container-level provenance recording. All nine mathematical/helper function ASTs are identical. `portable.diff` and `portable_preparation.json` record that comparison.

Testing used a relocated minimal archive subset with no `.git` directory and an empty process `PATH`. Both the default compressed invocation and explicit raw-result invocation passed. All scientific output fields equal the original observed review. The only excluded comparison fields were input path/container hashes, verifier hash, start time and runtime. Both invocations verified 26 radii, 72 row deletions and 220 stored search-event certificates. Optimized Python, an existing output, and wrong decoded result content were rejected before any output was written. The test used CPython 3.12.14 on Windows; other platforms have not been executed in this review.

| Artifact | SHA-256 |
| --- | --- |
| Original verifier, preserved | `e971da5fab874f0843b58d2e85dc9fc6b28939254eb4ee3c19d73466f46bc46f` |
| Portable verifier | `d7d608d207509b5ccd87394d9fc37131611683698226d29901746ece1e72fd1c` |
| Portable test receipt | `5cdc6707aa66e9e9b4fd2896976b47ef010bd5e0b8f160865a381b16660a7439` |

For archive integration retain the original and portable verifiers, the original review/receipt, this document, `portable.diff`, `portable_preparation.json` and `portable_receipt.json`. The fixture and duplicate replay-output directories are disposable review evidence, not required archive inputs. No manuscript or publishing asset is needed. The numerical and interpretation limits in `review.md` remain unchanged.
