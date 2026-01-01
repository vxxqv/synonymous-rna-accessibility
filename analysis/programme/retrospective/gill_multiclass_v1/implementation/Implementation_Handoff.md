# Frozen Gill multiclass implementation handoff

**Ready for root's independent implementation review. No real training or validation was executed and the observed workbook was not opened.** The source and tests are held unchanged at the hashes below. Root must review and commit them before permitting new training counts, then commit the resulting training receipt/predictions before validation.

## Portable files and verification

- `gill_multiclass.py`: portable train/validate CLI; NumPy and Python standard library only. All input paths are explicit; no Git checkout, network or author notebook is needed.
- `test_gill_multiclass.py`: synthetic-only tests, never opens the observed workbook.
- `metadata_preflight.json`: actual frozen protocol/inventory geometry check, with zero outcome exposure.

The final run passed **27 pytest cases in 1.70 seconds** on Python 3.12.14 / NumPy 2.5.4. Cases cover poisoned D headers and values in train, poisoned unselected rows, forbidden B/C/E/F values in a validation reader fixture, invalid counts, missing/duplicate selected genotypes, formula rejection, zero counts, all-zero residuals and undefined ranks/cosines, exact integer ranks above floating-point integer resolution, and five training-gate tampering cases that prove the D reader is not called first. A copied script completed train and validate in a temporary directory with no Git repository. Optimized Python is rejected by the CLI.

Independent numerical references include OLS via least squares for categorical projection, a connected nonclique codon factor, 84 random/tied/negative-score Floyd max-min comparisons against the component-merge maximum-barrier algorithm, and 31 direct quadratic-decomposition checks including a reused DNA denominator with zero RNA signal. The actual frozen sequence inventory passes membership, ranks, Cartesian-factor and graph-edge checks for 2,176 classes / 38,338 vertices, partitioned 1,124 / 1,052 classes and 20,028 / 18,310 vertices. This metadata check opens no source workbook.

## Access and freeze gates

`train` first verifies the externally supplied protocol hash, decompressed inventory hash, selected membership and source-workbook hash. Its reader interprets A first; only selected genotypes have B/C/E/F/G interpreted. It never interprets D, including its header. Unselected rows have only A interpreted. Formula, error, boolean and noninteger/nonfinite/negative selected counts fail. Zero RNA and individual zero DNA are retained. Row positions, complete fixed key membership and duplicate protection are enforced.

Training writes deterministic compressed `counts.json.gz` and `predictions.json.gz`, `training_summary.json` and a receipt containing source/protocol/inventory/code and artifact hashes. Predictions preserve the fixed fields, categorical projection, residuals, raw counts/logs and training barriers for all three pseudocounts. Raw-integer rank diagnostics use the actual count ordering, not rounded logarithms.

`validate` requires the training receipt SHA-256 as an explicit pinned argument. Before calling the D reader it verifies that pin, all input/code/source hashes, the exact artifact set and artifact hashes, saved count membership/types, and reproduction of the saved training computations. It then reads only A/D/G for the fixed keys and checks G against the historically saved counts. This is a local immutable-artifact gate. The required historical Git commits are verified externally by root; the portable script does not falsely claim it observed a remote commit.

Validation preserves the held-out counts, per-class metrics and all summaries, all three pseudocounts and both partitions. It reports the primary weighted residual cosine and threshold classification, fixed residual prediction losses and whole-field centered shape losses against zero/constant/codon baselines, raw RNA and DNA rank diagnostics, RNA-only and DNA-only residual pooled agreement, all nine scenarios and all 45 quadratic terms, maximum barriers and fitted-codon barriers separated by clique versus nonclique architecture, every leave-one-stratum-out cosine and the maximum stratum energy shares. Undefined values are JSON null and retained with counts; no class is dropped because of its result.

## Numerical details to review

The categorical projector uses uniform conditional means on the four complete codon factors, P f=sum_j E[f|codon_j]-3 E[f]. Residual entries no larger than `128 * machine_epsilon * n * max(1,maxabs(field))` are set to zero as a fixed roundoff rule. This is documented in output and is not a biological threshold. Linearity and quadratic identities have scale-aware numerical checks. The primary residual is assembled from the same projected RNA/DNA vectors used by the decomposition, avoiding inconsistent denominator definitions.

The maximum barrier routine activates vertices in descending score order and merges their active components. Each merge at level t joins two components whose maximal endpoint scores are a and b; it covers a candidate barrier min(a,b)-t. Every endpoint pair first becomes connected at some such merge, and choosing component peaks maximizes its possible endpoint-relative barrier. Consequently the maximum over merges equals the maximum over all endpoint pairs. Ties yield zero-height merges and are allowed. Reported pair counts are called `endpoint_pairs_covered`, since the algorithm proves coverage rather than individually enumerating each pair. The independent Floyd tests check this specialization on synthetic inputs.

Raw maximum barriers are retained. Each barrier now includes a numerical tolerance `128 * machine_epsilon * n * max(1,maxabs(field))`, a resolved-positive flag and an unresolved-tiny-positive flag. Aggregate raw-positive, resolved-positive and unresolved-tiny-positive counts are separate. This is a conservative precision diagnostic, not a statistical or measurement-error threshold; no observed maximum is silently zeroed. A 180-field additive-clique stress check passes, and an explicit tiny positive remains in the output with an unresolved flag. Whole-field losses separately center training and validation fields: these are shape comparisons, not out-of-sample level predictions. The primary residual endpoint and every raw scientific metric are unchanged by this reporting amendment, checked exactly on a synthetic before/after replay. No stochastic confidence interval, biological noise ceiling or independence assumption is introduced. The failed reference-protein pilot stays quarantined. A positive barrier in a nonclique codon factor is not automatically an interaction between codons. The protocol's scale and biological limitations are preserved.

## Reproduction and execution order

Run the synthetic checks from the workspace root:

`work/env/Scripts/python.exe -m pytest work/gill_multiclass_cycle3/implementation/test_gill_multiclass.py -q`

Portable CLI templates, after root review and implementation commit:

`python gill_multiclass.py train --protocol PROTOCOL.json --inventory INVENTORY.json.gz --workbook SOURCE.xlsx --expected-protocol-sha256 4ab4e85c14203c0ab8c983e8433c766f1566f373a335fd3fda4fb3981371c1aa --output FRESH_TRAIN_DIR`

After root fixes and commits that training receipt/predictions:

`python gill_multiclass.py validate --protocol PROTOCOL.json --inventory INVENTORY.json.gz --workbook SOURCE.xlsx --expected-protocol-sha256 4ab4e85c14203c0ab8c983e8433c766f1566f373a335fd3fda4fb3981371c1aa --training FRESH_TRAIN_DIR --expected-training-receipt-sha256 PINNED_RECEIPT_SHA256 --output FRESH_VALIDATION_DIR`

Each output directory must be new. The script does not overwrite earlier results or resolve a mismatch by choosing different classes or settings.

## Exact held hashes

- `gill_multiclass.py`: `9d60eabe28cf8a8f6c37f17a39c1adcbec049eca04e526c04751605657569268`.
- `test_gill_multiclass.py`: `4874bb513005084b36a1a10cfd548f3f3743f0e69c14297de4cc599bf42fd570`.
- `metadata_preflight.json`: `31c7622d71dd447a0ea410943126ce993e5b07fcfba45b7c58b2899e7e5bf054`.

Protocol commit supplied by root: `d640324987908a87fa6ce9ef60288d80e30da032`. No actual outcome cells, manuscript/repository edits, commits, publishing or physical interventions occurred in this implementation task.

## Preserved pre-outcome amendment history

The initial 24-test source SHA-256 was `3faca97fff160a7be0e255bdcc102cc79a8b0ab47d35f4e30d5af9f6ae04f146`; its test receipt SHA-256 is `59d255fcf9a881afa063071d3e592e73c8dcbcd3f87570e8e1fee0373e46779c`. Both source, tests, receipt and handoff remain under `history/pre_precision_reporting/`. The root requested the precision-count and shape-label amendment before outcomes. Its first synthetic run caught a NumPy boolean JSON-serialization issue, which was fixed by converting diagnostic flags to native Python booleans. The final 27-test run passes. `precision_amendment_parity.json` records 48 exactly equal scientific summary comparisons and exact equality of all raw per-class metrics after removing the added flags. No real outcome was accessed in either version.
