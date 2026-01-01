# Offline replay of the final CNOT3 numerical verification

This wrapper replays the existing independent verifier on the archived frozen inputs. It does not run the author pipeline, alter models, choose new penalties, or access an XLSX workbook. The unchanged independent verifier does recompute the fixed fold calculations and declared predictions, then checks the saved prediction and evaluation artifacts. This is numerical verification, not a new analysis or biological validation.

## Requirements

Use an already available environment with CPython 3.12.14, NumPy 2.5.4 and Biopython 1.88. The wrapper uses only the Python standard library; the unchanged numerical verifier requires NumPy and Biopython. It rejects a different runtime. No installation, Git executable, repository history, credentials, network connection or external source data are required.

Keep `replay.py` beside this `manifest.json`. The wrapper hard-pins the manifest SHA256 and checks all 24 listed files before any scientific JSON is decoded. The six cohort review files are included as integrity-checked provenance, not as a repeated cohort computation. All paths in the manifest are relative to the archive root and match the integrated archive structure.

Suggested placement inside the computational archive:

```text
analysis/programme/retrospective/cnot3_v2_exploratory/final_numerical_replay/
    replay.py
    manifest.json
    README.md
    test_replay.py
```

From the archive root, using the compatible interpreter:

```text
python -I -B analysis/programme/retrospective/cnot3_v2_exploratory/final_numerical_replay/replay.py --archive-root . --out ../cnot3-final-replay
```

The output path must be a new directory outside the archive. The command copies the two immutable verifier files byte for byte into that directory and executes them with explicit archive input paths, isolated Python mode and bytecode writing disabled. It never imports the author implementation. The author source is required only for its frozen hash.

Optional guard tests require only the standard library:

```text
python -I -B analysis/programme/retrospective/cnot3_v2_exploratory/final_numerical_replay/test_replay.py
```

## Expected result and acceptance

The final verification has 59,268 checks and 240,144 numerical comparisons, covering nine prediction vectors, 35 out-of-fold vectors and seven ridge-penalty selections, together with the final evaluation. On the tested Windows environment, replay produced a byte-identical `verification/verification.json` and the same maximum absolute discrepancy, `3.397282455352979e-13`.

Acceptance uses the unchanged verifier tolerances: prediction absolute `1e-9`, relative `1e-10`; score absolute `1e-12`, relative `1e-10`. The historical maximum discrepancy is diagnostic, not a tighter acceptance threshold. The wrapper requires the original counters, exact discrete fields, original tolerances, no failures and unchanged input hashes. It reports whether the complete receipt and its bytes match the archived receipt exactly. It does not silently permit a changed input or runtime.

Inspect `replay_receipt.json`, the unchanged verifier's `verification/verification.json`, and `stdout.txt`/`stderr.txt`. A failure preserves its output and `replay_failure.json`. Nothing is overwritten on retry; choose a fresh output directory.

## Scope and limits

This replay was tested in an isolated directory containing precisely the 24 pinned archive files and no `.git` directory. It made no Git or network call. This is code-level offline operation, not a claim of an operating-system network sandbox. Windows long paths are handled explicitly by the wrapper. Other operating systems and BLAS builds were not tested; the existing verifier's version and numerical guards remain in force.

The saved prediction commit receipt is checked as an artifact. Historical Git ordering is not independently reauthenticated without repository history or a remote service. Source eligibility, grouping, and cohort construction are not rerun here; their independent review evidence is pinned separately. The replay does not resolve the published cohort-size discrepancy, decay measurement calibration, experimental dependence, causality, or biological generalization.

## Integration boundary

Copy only the four files named above into the suggested directory. Preserve the existing verifier and all frozen inputs byte for byte. The complete manifest also requires `evaluation/evaluation.json`, `evaluation/receipt.json`, `independent_evaluation/verification.json`, `independent_evaluation_authorization.json`, `operational_prediction_commit_receipt.json`, and the existing cohort, prediction and independent-prediction artifacts. Cohort review evidence belongs under `analysis/programme/retrospective/cnot3_v2_exploratory/cohort_review/`.

The local preparation fixture, absolute-path copy map and earlier candidate manifest are not runtime inputs and should not be copied into the computational archive. The isolated replay receipt and this review report may be archived separately as verification evidence.
