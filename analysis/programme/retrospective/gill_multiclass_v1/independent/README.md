# Offline independent Gill verification replay

This replays an independently authored, pre-outcome frozen numerical verifier and its final complete comparison. It uses the archived repository layout and makes no network or Git calls. It does not open the source workbook, FASTQ, figures or manuscripts. All scientific inputs must already exist at their manifest-pinned repository paths. NumPy is the only third-party prerequisite; no dependencies are vendored or installed by the replay.

Run with **CPython3.12.14 and NumPy2.5.4**, from any directory:

Python `-O` and `-OO` modes are explicitly rejected by a `RuntimeError` before imports, scientific input reads or output creation, because the frozen scientific implementations require assertions. The original wrapper and its package manifest are preserved in `orchestration_history/`. This guard is the only execution-logic change after the successful isolated full replay; its refusal behavior is tested separately without recomputing the mathematics. The full replay receipt therefore records the preceding wrapper revision, while the guard-check receipt records the current revision.

```text
python /path/to/repo/analysis/programme/retrospective/gill_multiclass_v1/independent/replay.py --repo-root /path/to/repo --output /path/to/new-independent-replay
```

On Windows, quoted absolute paths work with the same arguments. `--output` must not exist and must be outside the archived input repository. The runner refuses an existing output, authenticates all source/history files and external inputs, computes the independent results before author comparison, checks all3,352,385 comparisons, and checks that inputs stayed unchanged. A failed check stops execution and retains partial evidence. No overwriting, input rewriting, selection or tolerance fitting is performed.

The repository inputs are referenced, not duplicated here: `analysis/programme/retrospective/gill_multiclass_v1/{protocol.json,inventory.json.gz}` and `results/programme/gill_multiclass_v1/{train,validate}/` files enumerated in `package_manifest.json`. The package works from a normal extracted archive containing that layout. It does not depend on the original user's workspace path. Results and historical records are kept separately under `results/programme/gill_multiclass_v1/independent/` in the delivery tree.

The byte-preserved scientific evaluator, comparator, design and freeze are in `frozen/`. The evaluator's original directory globals are rebound through a small orchestration adapter. The comparator is executed with only explicit path bindings and the current replay output-integrity digest substituted in memory. The exact10 archived path tokens map `training_v1/` to `train/` and `validation_v1/` to `validate/`. Its formulas, membership, pseudocounts, graph coverage, budgets and discrepancy logic are unchanged. `binding_record.json` records the executed binding digest and historical versus replay result digests. The frozen scientific source itself is never rewritten.

Float64 results can differ by tiny rounding across supported numerical backends. The original independent result's byte digest is retained and byte identity is reported separately. A fresh result is authenticated by its current digest, then every scientific quantity is checked against the fixed author arrays using the pre-outcome numerical budgets. Expected coverage and all scientific headlines are also checked. There is no claim that numerical backend changes preserve compressed result bytes.

The first comparator's representation mistakes remain in the history: raw Gram products were initially compared with class means, and named library pairs with stored index pairs. The first failed comparison JSON is losslessly gzip-compressed; its original uncompressed digest is recorded. Both prior adapter snapshots are byte-preserved. These schema corrections changed no mathematical formulas, membership or tolerances. The earlier failed metadata-stratum assumption is also retained with its exact failed source. The recovered phase1 review matches its freeze-pinned original digest.

Successful replay means agreement on processed exports, not biological validation. The fixed0.5 target fails in both partitions at all3 pseudocounts. Primary validation cosine is0.07787615736918292, and fixed residual prediction has negative gain against zero. Local classes share one truncated reporter. Neither this replay nor its coverage establishes biological replicate independence, calibrated noise, causal interactions or independent-gene generalization.

**Historical Git and access order are not reauthenticated.** Original commit IDs and access-order statements remain provenance assertions from the original receipts. The offline replay authenticates archived file bytes and recalculates their numerical consequences. It does not independently reconstruct who accessed data when or query a Git history.
