# Portable replay of the frozen observation-cost audit

This package replays the already reviewed 384-case exact finite study. It introduces no new graph, policy, prior, price, mathematical method or outcome. The original source, protocol, freeze, results, interpretation and independent review are preserved byte-for-byte in `frozen_sources.zip`. The large `cases.json` is losslessly compressed inside that archive; extraction restores its original bytes and SHA256.

## Prerequisites and command

Use Python 3.12 with **NetworkX 3.7 already installed** in that Python environment. These are the versions used for the successful replay: CPython 3.12.14 and NetworkX 3.7. `requirements.txt` pins the only external package. Dependencies are not bundled or installed by replay. No Git, network access, original workspace, repository checkout, manuscript or biological data is needed.

Keep these six distribution files together: `replay.py`, `README.md`, `requirements.txt`, `correlated_prior_addendum.md`, `frozen_sources.zip` and `integration_manifest.json`. After extraction, run:

```text
python -I -B replay.py --out NEW_DIRECTORY
```

`NEW_DIRECTORY` must not exist. The isolated-mode flag `-I` ignores Python environment-path customization; `-B` disables bytecode writing. The entrypoint also disables bytecode itself and rejects `-O` or `-OO`, because the frozen checks use assertions. It checks the installed NetworkX version and makes no installation or network call. A replay normally takes tens of seconds. The original author loop retains its 180-second cap and the separate reference process has a 240-second timeout. Failures return a nonzero exit code and retain an error receipt rather than changing results or selecting another case.

`receipt.json` in the new directory records exact comparisons, software, hardware, elapsed times, output hashes, and immutable input hashes before and after. It explicitly records `historical_git_verified: false`. This reproduces source and calculations without claiming to reauthenticate the historical Git events recorded in preserved manifests.

## Preserved sources and separate outputs

`integration_manifest.json` lists every archived original file and its hash. The imported selector `src/acquisition.py` remains SHA256 `6005d30d709dc8100a1de64822d8e80cf8a8dbc8d9c3a18de8d21f69622be6ee`; the author source remains `a82e32ab46bdc941235623b80ced7f600078f74f280d9aaeea0baf5d9aa2f107`. The original freeze remains `3686acf9a6f5c1210cff47630de0959f536d1c552d57577f7157c939bae8a9d3`.

The new output separates:

- `preserved/`: immutable original files, including all historical receipts and the original large case file.
- `reconstruction/`: only static source and contract inputs are copied here before any run. Fresh author and reviewer results are created in their expected relative layout. Historical output files are never placed in active output paths before execution.

The adapter invokes the unchanged author model, Bellman solver, world-policy evaluator, C4 exhaustive-tree check and summary functions. Its outer case loop preserves the original order, bounds, checks, runtime cap and serialization. It avoids the original `run()` function's platform-specific relative-path handling and Git-dependent manifest construction. Fresh orchestration provenance is in the portable receipt rather than a forged historical manifest.

The unchanged, separately authored reviewer executes in a separate Python process. It uses integer world-support masses, its own graph traversal and exhaustive subset separators, with no call to the creator solver, acquisition module or NetworkX. The only compatibility adapter intercepts parsing of the exact frozen freeze text and normalizes its relative path separators in memory. No other JSON is altered, and the freeze file and reviewer source remain byte-identical. The reviewer still checks the hard-coded freeze hash and all frozen source hashes. Its own input hashes are checked against actual files inside the reconstruction. The review's fixed-grid summaries are also reaggregated from its independently recomputed per-case expectations.

## Exact comparisons and provenance exclusions

Every field in fresh `cases.json`, `summary.json` and `c4_bruteforce.json` is compared exactly with its preserved counterpart. The decompressed bytes of `bellman_certificates.json.gz` must agree exactly. Raw gzip-header/compression bytes are reported separately because library and host serialization details may differ; original compressed bytes are always preserved.

The reviewer result is compared exactly except for these explicitly identified top-level provenance fields: `inputs_sha256`, `command`, `python`, `runtime_seconds` and `started_at`. Absolute paths and host text serialization can change those hashes; the fresh hash entries are independently verified against actual reconstructed files. Every count, cost, action choice, boundary interpretation, limitation, independence statement, status and per-case comparison must agree. No scientific field is excluded. The separately saved `summary_verification.json` must agree in full.

Historical manifests may contain old machine paths as provenance, but no replay dependency follows those paths. No source is rewritten, no original output is overwritten, and no Git query occurs. No PDF, figure, publication file, biological measurement or bytecode cache is included.

## Mathematical scope and prose qualification

Read the preserved `prior_and_proof.md`, `interpretation.md` and `reviewer/review.md`, together with `correlated_prior_addendum.md`. The addendum clarifies the generic correlated-prior paragraph: conditional expectations are defined at positive-probability histories, while universal stopping certificates must cover **all logically possible compatible worlds**, including zero-prior worlds. A zero-mass branch can be omitted from an expectation but still needs a sound counterfactual continuation if correctness is required in every world.

The experiment itself has full support at all three frozen probabilities, so that clarification changes no calculation. Its 384 cases are exactly 64 price vectors, three independent Bernoulli probabilities and two action menus. The 98,304 full-world policy trajectories and all 279,936 state values are finite exact evidence. Optional batches help in 184 of 192 matched cases and tie in eight; failures and all other comparisons are retained. The simple price-aware policy remains suboptimal in every case. These are summaries of this chosen finite grid, not population estimates, assay economics, noisy-observation results, algorithmic novelty or scalable optimization claims.

## Integration

Archive the six distribution files together. A successful portable `receipt.json` may accompany them as execution evidence. The duplicate extracted trees need not be archived again. The separate addendum must remain with the historical proof so its zero-probability qualification is not lost. Existing failed acquisition or Gill studies are neither changed nor reinterpreted by this package.
