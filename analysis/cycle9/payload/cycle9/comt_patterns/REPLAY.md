# Local comparison and direct reconstruction

Standard-library CPython3.12.14 was used. No Git/network access, R packages, new source extraction, fitting or source-control values are required for the comparison. The original source reconstruction remains separately conditional on its runtime and source contract.

Preserve the existing `observed_v1`, `reconstruction_v1` and synthetic receipt. A fresh reviewer can check the saved result with explicit paths:

```powershell
work/env/Scripts/python.exe -B work/research_cycle9/comt_patterns/verify.py --input work/research_cycle8/comt_results/actual_v3/normalize/contrasts.json --result work/research_cycle9/comt_patterns/observed_v1/result.json --inventory work/research_cycle9/comt_patterns/metadata_inventory.json --freeze work/research_cycle9/comt_patterns/execution_freeze.json --out NEW_REVIEW_DIRECTORY
```

The checker imports no production analysis function. It uses separability to reconstruct the ranges, direct enumeration to verify every joint record, and sets of exact signs to classify them. The same-author authorship is disclosed; a fresh reviewer must judge the input and implementation contract independently.

For an explicitly authorized clean-copy replay, retain the files named in `execution_freeze.json` together, and supply byte-identical input, inventory and contract paths. The implementation refuses overwrite, changed code/input hashes and optimized Python before outcome decoding:

```powershell
python -B analyze.py --freeze execution_freeze.json --input PATH_TO_SEALED_CONTRASTS_JSON --inventory metadata_inventory.json --contract PATH_TO_SELECTIVE_CONTRACT_JSON --out NEW_RESULT_DIRECTORY
```

The second command documents reproducibility; it was executed once on the observed source in this task. Do not erase or relabel the original freeze/result as a new prospective validation. The comparative outcome was already exposed in Cycle8 before this exploratory follow-up was proposed.

The compact CSVs are presentation exports of the complete pinned result. Exact fractions remain authoritative; CSV decimals are displays. `metadata_inventory.json` retains all53 sites and all125 IDs, including singleton sites. Synthetic receipts are historical and immutable; rerun synthetic checks only in a disposable copy because their output also has an existence guard.
