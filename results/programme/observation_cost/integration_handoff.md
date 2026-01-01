# Observation-cost portable replay handoff

Status: **passed, finite implementation/result reproduction in the stated range**. This packaging task added no mathematical method, graph, policy, experiment, parameter or case selection. No original study, repository or manuscript file was edited.

`Observation_Cost_Portable_Replay.zip` is the ready distribution. It contains exactly six files: `replay.py`, `README.md`, `requirements.txt`, `correlated_prior_addendum.md`, `frozen_sources.zip` and `integration_manifest.json`. Its loose copies are identical. It requires preinstalled NetworkX 3.7 and Python 3.12; no dependency is vendored, installed or retrieved by replay. The command is `python -I -B replay.py --out NEW_DIRECTORY`.

## Verified replay

Only those six files were copied into `qa/no_git/package/`. A single full replay ran into the previously nonexistent `qa/no_git/replay/` directory. No Git metadata exists in that test directory or its ancestors. Its `receipt.json` is SHA256 `ce74ee6562e04f6e5ddb311ddbd6db5351832633f75945dd6400de1e5528d2e1`. It records `status: passed`, no network or Git use, and `historical_git_verified: false`.

The creator calculation took 7.35 seconds; the separate reference process took 7.64 seconds. CPython 3.12.14 and NetworkX 3.7 were observed. The following scientific outputs agreed exactly:

- Full JSON contents of all 384 cases, the summary and 24 C4 exhaustive-tree checks.
- Decompressed exact bytes of all 384 Bellman certificate tables. Raw gzip bytes also agreed on this host; only wrapper/header compression metadata is allowed to differ on another host.
- Every scientific field of the separately authored review, including 279,936 state values, 134,784 selected actions, 1,536 expectations, 98,304 world trajectories, 36,868 full batch charges and 10,384 partially observed batch actions.
- Independent exact reaggregation of every scientific summary field, plus full equality with the preserved `reviewer/summary_verification.json`.

The creator additionally reproduced all 492,480 Bellman action inequalities, 4,096 compatible state/completion certificate checks and all 729 partial-state labels. There was no partial grid, mismatch, case exclusion or result adjustment.

The review receipt excludes only top-level `inputs_sha256`, `command`, `python`, `runtime_seconds` and `started_at` from the historical/fresh comparison. These are provenance, not costs, choices or other scientific values. Fresh reviewer input hashes were independently checked against actual reconstructed files. No other scientific field is ignored.

All 16 original payload files and all eight static active inputs remained byte-identical. All four freeze input hashes matched. `qa/validation.json`, SHA256 `82790bad137a83697fd9560a4017226b6c9758965c4838b1a7dc5b9b20b430f6`, additionally records successful refusal of an existing output directory and disabled assertions, no forbidden asset or bytecode, identical copied distribution files, and unchanged original sources/results. `qa/check_package.py` performs those bounded checks without rerunning mathematics.

## Compression and integration hashes

The internal source archive is 753,462 bytes. The original `cases.json` is 18,565,256 bytes, compressed losslessly to 381,535 bytes within the ZIP. Its restored SHA256 remains `bfc3259d7748a18532033fecd3545e9f56849fa8c07cd666a9ed17549b5f3d9c`.

| Distribution file | SHA256 |
|---|---|
| `replay.py` | `bf4c5ddd53f9f7e8580554f4b0f4e27de67019009e1daa940399b79caf160985` |
| `README.md` | `629bd5712d4108988a6423ffd58ccd5b24a6af04c545c1c046ca902f97092df9` |
| `requirements.txt` | `63b83d49071ae093c187423c768c6711a3716cda354d3a24e2bb10b0f536d34c` |
| `correlated_prior_addendum.md` | `6aab759cf05cfea15b5fae336e3dd0a811aa570bcb0f1266cc3a3fb87cec527c` |
| `frozen_sources.zip` | `7db60195b976b22228af281d38277c0c13421be2471749772fdce0824b3a2e70` |
| `integration_manifest.json` | `4b9c3439168fea7c3fa9603b24668211e5132c426875e87b7ef5053f41571c6c` |

The preserved creator source is `a82e32ab46bdc941235623b80ced7f600078f74f280d9aaeea0baf5d9aa2f107`; acquisition source is `6005d30d709dc8100a1de64822d8e80cf8a8dbc8d9c3a18de8d21f69622be6ee`; independent reference is `30dcc83f2a5973211baee0d4c41c51c2c4779154c8d8f1466c5ad12062bc809c`. Every preserved file is individually hashed in the integration manifest.

## Orchestration and proof boundary

The original relative layout is reconstructed. Historical results stay only in `preserved/`; static inputs alone are copied before active runs. Original mathematical functions, case ordering, action conventions, checks, serialization and the author runtime cap are unchanged. The portable wrapper replaces historical platform-specific path and Git manifest orchestration with explicit fresh provenance.

The unchanged independent reviewer executes separately. Its only compatibility adapter normalizes relative freeze-path separators in memory when parsing the exact frozen freeze text. It does not rewrite that file, reviewer source or any scientific JSON. Old absolute commands and paths remain inside historical receipts as provenance, never as active dependencies. NetworkX is a pinned preinstalled prerequisite rather than a newly packaged dependency tree.

The separate prose addendum is mandatory for integration with the proof. It states that conditional expected costs under a generic correlated prior are defined on positive-mass histories. Universal stopping certificates still cover all logically compatible worlds, including zero-prior worlds; omitting a zero-mass branch from an expectation does not certify its answer. A sound counterfactual continuation can observe remaining singleton coordinates. The frozen experiment is full-support and numerically unchanged. No correlated-prior computation or new uncertainty claim was added.

No primary PDF, publication export, manuscript, figure, biological measurement, bytecode cache or Git metadata is distributed. The six distribution files and the successful receipt suffice for integration; the large generated duplicate extraction trees do not need to be archived again. Root retains responsibility for review, integration, commits and publication. This task stops after the single passing full replay and bounded integrity checks.
