# A3 portable replay handoff

Status: passed. The archive-ready distribution is `A3_Portable_Replay.zip`, containing only `replay.py`, `README.md`, `computational_sources.zip` and `integration_manifest.json`. The loose copies are identical. Run `python -B replay.py --out NEW_DIRECTORY` after extraction. Python's standard library suffices; CPython 3.12.14 was tested. No Git or network is used. No original source, repository or manuscript was edited.

## Successful isolated run

The four distribution files were copied alone into `qa/no_git/package/`, with no Git metadata in that directory or any ancestor. Its successful fresh output is `qa/no_git/replay/`. The receipt is `qa/no_git/replay/receipt.json`, SHA256 `fd6e50cb5adc4a273fe1b19eab9b503d88a87644278c3b40dc6da12f8fe17f32`. All three scientific comparisons passed exactly:

- Original author finite receipt: all scientific fields, counts, status and original input hashes agree; only runtime and software provenance are excluded. Its earlier full replay core also agrees exactly. The original 2,818 parameter inputs and 75 small MaxCut graph checks were replayed in 3.54 seconds.
- Separately authored reference: complete JSON equality, with no excluded fields. This reproduces 76,544 projection checks, nine additional one-factor checks and 217 scaled switch checks in 1.46 seconds. The 1,296 Q2 parameter inputs have 929 distinct fields; the 32 deterministic Q3 inputs have seven distinct fields.
- Corrected revision-2 suite: every mathematical, storage, order, bookkeeping, domain-restriction and coverage field agrees; only runtime and provenance are excluded. Runtime was 0.07 seconds. No historical Git query was rerun.

All three freezes were verified against preserved files, with 5, 14 and 17 input entries respectively. All 36 original payload files were byte-identical before and after packaging and replay. Static active copies were also unchanged. `qa/validation.json`, SHA256 `f4812d083b3a6a0b5435ad57fc6594a95a3cc91d9ad70120828c0bca23ee6004`, records those checks plus successful rejection of optimized Python and an existing output directory. No bytecode cache was produced. `qa/check_package.py` is the bounded integrity/refusal checker; it does not rerun mathematics.

## Integration hashes

| Distribution file | SHA256 |
|---|---|
| `replay.py` | `3acdb492d2d08c298faf071ab26b69549312f95eaa83611ade8aa43d3bd38911` |
| `README.md` | `3846c46b080baf49e1997f790f97b1aa21c667ac25043efcc0483a8717487cbe` |
| `computational_sources.zip` | `d62516c7a665971408e7ebd37c8a1f9e07170f29ed37ef134c9c0f89b52527e5` |
| `integration_manifest.json` | `d5ddf25139ddab48c56b9a82cf17296aaab8fae6ec53ed21a7d84242ddf49181` |

The manifest enumerates every original source, receipt and correction file. The creator source remains SHA256 `1ceeb89a3c40c2bae2c4995c8014c6163f51e8eff0be75de5c48ee988c5179e2`, and the additive dependency remains SHA256 `96c6d35a6243dcf1ba0dbd08100d82dbbf553a583e3844e9ac6a31646d5ec63a`.

## Faithful scope and correction history

The current proof is `revision_2/claim_and_proof.md`; the current implementation cost statement is `revision_2/complexity_contract.md`. The correction ledger and trace erratum are mandatory companion records. Historical proof and review text retain the superseded total-storage reading as history. The initial diagnostic retains its two inflated trace counters as history. `results_v2.json` is the corrected accounting receipt. No obsolete complexity claim is promoted to the current contract.

The portable adapter imports unchanged author functions in their original invocation order and imports the unchanged corrected suite. It normalizes path separators only when interpreting historical freeze keys. It separately executes the original independent reference. It bypasses old platform-specific path orchestration and Git-dependent receipt construction, not mathematical code. Historical manifests may contain old machine paths as preserved provenance; no replay dependency follows those paths. The fresh receipt explicitly records `historical_git_verified: false`, preserving the distinction between source-hash reproduction and reauthentication of historical repository events.

The archive contains no literature PDF, publication export, manuscript, figure, expression data, Git tree or bytecode cache. The small historical PDF-retrieval receipt remains metadata only. Finite agreement is corroboration, not a universal proof or novelty demonstration. Cubic-MaxCut hardness retains its established external-theorem dependency. Baseline-component trimming remains caller preprocessing, natural-order width remains the implemented convention, and full retained witnesses and envelopes remain in the corrected total-storage bound. No new mathematical method or biological claim was introduced.

Only the four distribution files need integration. The successful receipt and this handoff may accompany them; the generated duplicate trees under `qa/no_git/replay/` need not be included in the computational archive. Root retains responsibility for integration, commit and release. This task stops after the passing portable replay.
