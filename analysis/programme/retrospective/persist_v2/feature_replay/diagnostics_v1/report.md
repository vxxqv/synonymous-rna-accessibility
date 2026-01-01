# Linux feature replay: failure visibility, not a diagnosed scientific error

The child cause is **not recoverable from the retained evidence inspected**. GitHub run40178313369, job121123714821, exposes only the parent exception `Worker stdout/stderr is nonempty; stop without retry: 310030T3`. The run has no retained artifacts. No actual worker stderr, stdout or exit code appears in the downloaded job log. This report does not identify a cause or claim that the first fold was reached.

Repository review point: `602ac0559efcc9c5e10410d85f6a164761dc2b3f`. The decoded job log is archived as `job_121123714821.log`; the independently retrieved artifact listing is `github_artifacts.json`. The earlier Windows full-feature receipt remains unchanged: all37 rows,3071 feature values,3367 floating comparisons, zero maximum absolute error and zero warning rows. Linux failure is a separate unresolved portability gate. Steps preceding feature replay passed in this Linux job; the five later scientific checks were skipped after its failure.

## Why the exception is hidden

In archived `analysis/programme/retrospective/persist_v2/feature_replay/replay.py`, line105 writes each subprocess's decoded streams to `<output>/<SequenceID>/stdout.txt` and `stderr.txt`. Line106 then calls `accept_process`. Lines31–34 reject nonempty streams before checking the child return code, without printing their content. A nonempty stream can accompany either a failing or a successful process, so the visible exception does not establish a child exception or its exit status. Line126 saves a parent `FAILURE.json`, including completed IDs and the parent traceback.

Current `.github/workflows/verify.yml`, lines40–41, runs the replay into `${{ runner.temp }}/ci-persist-features`. It neither prints these generated stream files nor uploads their directory after failure. The missing evidence was on the ephemeral runner filesystem. Re-downloading the same retained job log cannot recreate it.

## Minimal recovery without changing scientific behavior

`workflow_diagnostic_fragment.yml` proposes only an ID on the existing feature step and two subsequent failure-conditioned observation steps. `observe_failure.py` prints the already-written parent metadata and worker stream/warning diagnostics as JSON-escaped text. It records SHA256 and length, including empty streams, and only hashes parameter exports. JSON escaping keeps any embedded GitHub workflow-command strings inert. The final upload preserves the entire generated feature-output directory, including complete stored stream files. The upload runs with `always()` if the feature step failed, even if printing also fails.

The original failing step remains failed. There is no `continue-on-error`, warning suppression, scientific acceptance change, threshold change, worker retry or fallback. No frozen source file needs replacement. Integration of these snippets and any future CI execution are reserved to root. The snippets themselves do not recover the lost files from run40178313369; they make a future authorized run inspectable.

Important byte distinction: the frozen parent already decodes subprocess output as UTF-8 with replacement before writing its files. An upload preserves these stored diagnostic files exactly, not any undecodable original OS-pipe bytes that the parent has already replaced. The observer adds no further truncation. Its SHA256 values authenticate the stored files, not raw child pipe bytes.

## Optional first diagnostic with zero folds

`config_preflight.py` and `workflow_config_only_fragment.yml` provide a separate narrower probe if root prefers no new fold until configuration has been checked. The script pins the original protocol, author source and Turner2004 file hashes; requires the same CPython3.12, NumPy2.5.4, ViennaRNA2.7.2 and double-precision PF build; imports the frozen author; and calls only `configured_md`, `snapshot_md` and parameter export. Those are the same explicit pre-fold operations in replay lines68–72. It never calls `fold_one`, `sequence_features`, `fold_compound`, MFE or PF, reads no sequence metadata, feature matrix or response table, and cannot produce a feature-replay success receipt. Importing the pinned author is side-effect-free with respect to source tables and folds, as verified by source inspection.

The preflight records the failing stage, actual exception/traceback, model snapshot if reached, exported parameter hash if reached and before/after source hashes. Any Python warning is fatal. Output must be fresh and outside the repository. A successful preflight verifies configuration only. It neither explains the old failure nor proves the full Linux numerical replay. A failure is new actual configuration evidence, not retrospective proof of the lost worker's exact exception.

Proposed invocation after root integrates the script at the path in its snippet:

```text
python -I -B analysis/programme/retrospective/persist_v2/feature_replay/config_preflight.py --repo-root . --output /path/to/new-config-only-output
```

## Bounded validation and limits

No Linux runtime is already available locally: `wsl --list --quiet` reports WSL is not installed; Docker is not available on PATH. Nothing was installed. No real preflight, fold, feature reconstruction, model fit or outcome read was executed for this task.

Seven synthetic-only tests passed. They cover full stream/empty stream preservation at the stored-file level, JSON-safe log output, a strict diagnostic filename whitelist, explicit absent-output status, a fake configuration success with fold/feature functions that raise if called, actual fake exception retention, fatal synthetic warnings, model-snapshot and parameter-byte mismatch, output overwrite refusal, and early `-O`/`-OO` refusal for both scripts. They establish diagnostic behavior only. They do not validate Linux configuration or numerical reproducibility.

The shared scientific files and CI were read only. The owned `manifest.json` records inspected input hashes and all delivered diagnostic file hashes. Root's unrelated concurrent work is outside this audit. There is no confirmed cause and no proposed scientific correction.
