# Independent observed branch D verification

Verdict: **PASS. All2556 complete pair records, all69012 scenario capacities and all exact interval fields agree. No mismatch or incomplete scenario was found.**

Authorized primary SHA256: `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466`.

Authorized branch D result SHA256: `0c2ad1a6a602f80e44ff05eaf3321522d9eaf0940382c599309a04d33f4b383e`.

Frozen manifest SHA256: `340624e795e600abbbaeb3daefa89b054ceca7e356855014d320df20b043a10f`. A separate read-only Git-content check confirms this manifest and all six frozen source files exactly match revision `7eae1381c781c93e324f01d4e0d07692bceabea5`. This is recorded separately in source_commit_check.json; numerical replay itself needs no Git and does not claim to authenticate historical access order.

## Exact reconstruction

The separately authored verifier imports no production scenario/ordinal function. It reconstructs204 Hamming-one edges on all72 declared sequences, retains all27 complete positive scenario fields, and computes capacities through BFS at every distinct score threshold. It independently recomputes the vertexwise minimum and maximum fields, rather than substituting selected scenarios. The result is cross-checked with the previously independently verified Cycle6 primary and all27 original scenario capacity tables.

Every full record is compared, including scenario order, E, C_fix, C_ad, C_max, C_upper, exact/bounded status, all five interval records, lower/upper closure flags, empty intervals, exact widths and both obstruction union intervals. All10412 critical boundary and between-boundary membership checks agree. No endpoint/scenario/unknown exclusion or favorable pruning was used. The source has no unknown cells; an unexpected unknown fails this complete-source verification rather than being silently dropped.

| Quantity | Independently reconstructed count |
|---|---:|
| Fixed-route interval nonempty | 2556 |
| Adaptive-routes-only interval nonempty | 72 |
| Mixed-scenario-connectivity interval nonempty | 22 |
| All joint scenarios blocked while optimistic box connects | 0 |
| Upper-box obstruction interval nonempty | 0 |
| All-scenario obstruction union interval nonempty | 0 |
| Adaptive-only and mixed overlap | 13 |
| Adaptive-only or mixed union | 81 |
| C_fix equals E | 2475 |
| C_fix below E | 81 |

Pair counts above are not a categorical partition. Different thresholds for the same pair can fall in different categories. All pairs and all27 capacities appear in observed_run/all_pairs.csv, not only exceptional pairs. Its SHA256 is `62e4fe15fee44fd1b4d6fad086a32633fbb06111eea33c4680c5bb77a8fa50c5`.

## The2475-pair statement and tau0

For any pair, E is the largest common absolute threshold at which both endpoints are accepted in every scenario. A fixed path across all scenarios exists at threshold tau iff tau<=C_fix. Therefore C_fix=E is equivalent to existence of a path accepted across all27 fields at E, and that same path also works at every smaller eligible threshold. Thus the2475-pair statement is correct for the entire common-eligible interval[0,E], including its upper boundary. It is not merely a statement at tau0.

All2556 pairs have some fixed-route interval and accept a fixed route at tau0. This alone is weak: a connected graph with nonnegative scores always retains every vertex at0. Here every source score is strictly positive, so all2556 fixed-route intervals also have positive width; sufficiently low positive thresholds are likewise uninformative about high-score constraint. The2475 endpoint-limited criterion is the stronger exact statement. The other81 pairs fail to have a fixed route at E and equal the union of the72 adaptive-only and22 mixed pairs, with13 in both.

Zero all-scenario obstruction means there is no common-eligible threshold where every one of the27 joint fields is blocked. It does not mean every scenario connects: the mixed intervals explicitly show some do not. It also does not erase the earlier Cycle6 barriers computed at scenario-specific endpoint-relative thresholds, which use different eligible threshold domains.

## Quantifier proof and strongest alternatives

For each field r, C_r=max over paths of the smallest score on the path. Endpoints are included. Fixed route capacity equals max_P min_r min_v f_r(v), hence capacity of the vertexwise minimum field. Adaptive existence in every field has cut min_r C_r. Some field connects up to max_r C_r. The upper mixed-field capacity dominates every joint-field capacity. These give the ordered cuts and the stated partition, with weak >= acceptance assigning equality to the preceding interval.

The preflight's parallel-route diamond remains the strongest elementary alternative to a mechanistic interpretation: different fields can favor opposite intermediates while each independently admits a route. The serial-path example shows that a vertexwise optimistic box may synthesize a route absent from every actual joint field. Neither construction requires molecular interactions or a causal recovery effect. Actual field differences share samples/normalizations and are not independent biological confirmations. No practical absolute acceptance threshold was calibrated by this analysis.

The observed conclusions concern only27 supplied relative included-RNA enrichment fields, one fixed72-state published design, and common absolute thresholds. They do not establish true biological endpoint coverage by the empirical box, a molecular mechanism, measured evolutionary movement, function improvement, robust behavior under unmeasured bias, statistical significance or connectivity of the whole synonymous sequence space.

## Reproducibility and tests

verify_observed.py SHA256: `49445bb4fa6d22c268844c26a65169e4a2ef9480de0e6c05765f005eee891fd7`. It uses the standard library only, exact Fractions and no network/Git/production imports. Synthetic controls were run before observed execution: parallel/serial quantifiers, zero threshold, disconnected graph, ties, empty cells, a zero-width singleton and -O/-OO rejection.

A fresh archive-like layout without.git reproduced all scientific verification fields exactly and all_pairs.csv byte-for-byte using default repository-relative gzip inputs. Stored compressed-file hashes differ from the plain original input paths as expected, while decompressed hashes match exactly. The portable verifier always reports historical_git_verified=false. Guards reject existing output, output inside the checkout, missing root and an altered freeze before result output. No files in the original repository or authorized input files were changed.

The portable test harness initially pointed one directory too high when constructing its fixture source path; it failed before finding or reading an input. Correcting that test-harness path allowed the replay to pass. This was not a scientific verifier, numerical-result or threshold change. The standalone verifier bytes were fixed before synthetic/observed execution and remain unchanged.

The test harness portable_checks.py is workspace QA scaffolding and is not required for archival replay. The actual verifier accepts explicit paths or repository-relative defaults and has no workspace-layout dependency. README.md gives the exact archive/CI command. No new source/count outcome, broader branch, external action or commit was performed in this review.
