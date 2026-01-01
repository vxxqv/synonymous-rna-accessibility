# Bounded portable supplement adapter review

Verdict: **PASS**, with no blocking correction. Reviewed adapter SHA256 is `999f4b4ea8fc453daecfd7aa130412c1d321277cc298a85722e8cf4b6bed7aab` at `analysis/cycle5/replay_supplement_reviews.py` in the computational repository.

The adapter does not implement or change an arithmetic method. Lines 15-26 pin each separately authored oracle and its historical result, including the corrected reserved source. Lines 48-49 check these bytes before loading; lines 51-57 load that source and replace only its `REPO` and `HERE` locations. Both oracles import standard-library modules only, have no author-source import, and use `WORKSPACE` only to initialize the subsequently replaced `REPO`. Their existing input pins, methods and tolerances remain untouched. No new target or raw source is introduced.

Lines 58-65 remove only the `runtime` field, which the actual pinned oracles define as `platform.python_version()`, demand equality of all remaining parsed result values, then recheck the frozen source and historical result hashes. This is structural JSON equality, not a blanket claim that JSON serialization must be byte-identical on every runtime. Both existing full replay receipts in this runtime do in fact record result hashes identical to the historical result bytes. The associated saved replay outputs were checked directly without rerunning the scientific mains.

## Guard evidence

`check_adapter.py` passed 31 checks, retained in `checks.json`. Each arm was staged under a fresh temporary archive-like layout without Git or scientific input tables. Real command-line tests rejected altered oracle bytes, altered historical result bytes, existing outputs, an output within the checkout and the checkout itself. Existing output sentinels remained unchanged. Missing, file-valued and empty roots failed without creating output. Invalid arms failed. `-O` and `-OO` failed with the optimization guard before even resolving an intentionally missing input root.

A deliberately stubbed loader separately exercised location replacement and result comparison without executing either scientific main. Exact results and a changed Python-version field passed. A changed scientific field failed without a pass receipt. Source or historical-result drift during execution also failed without a pass receipt. These are adapter-control tests, not a new numerical reconstruction or independent rerun of the historical calculations. The final hash sweep confirmed every original reviewed file unchanged.

The test harness initially encountered Windows path-length limits in its own deeply nested temporary staging and cleanup. Only the owned harness was corrected to use extended-length staging and validated cleanup within its owned directory. The final 31-check run passed. This did not expose an adapter arithmetic defect or change any frozen source or result.

## Exact reviewed pins

| Arm | Oracle SHA256 | Historical result SHA256 |
|---|---|---|
| Development | `4db97846dba4e165a6cb63861ef369492b74dd531c726f8ef4cc65eee72cfd30` | `872b0ba211762aa1753117569126861fc4f32b19c1b377b554b3449d8a489cc4` |
| Reserved validation, correction v2 | `43cf8b1412bc571877dc8d63aa199b948e4df458f4a171977f4f56d80ac28311` | `7c231a7682ddc04c185bdaf3f9633cb0d6de83606b93092cdf4792c66d1a8bd3` |

Test source SHA256: `1ac6b583127258199517803f69455f3b59edb5274234a9926b6c429e1007c3b6`.

Check receipt SHA256: `00b8c4e51e545e1933fbfba17e9182899406dc288f32c4a86ec3318feb1446d7`.

No shared source edits, new outcome access, raw archive access, external action, commit or publication occurred. This review does not authenticate historical access order, repeat full graph geometry checks, or alter the failed biological validation interpretation. Root's existing two complete replays provide the scientific replay evidence; this bounded review checks the adapter and its meaningful failure paths.
