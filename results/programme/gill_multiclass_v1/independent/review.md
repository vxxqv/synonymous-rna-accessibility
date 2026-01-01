# Gill multiclass independent result review: prepared, validation unopened

**Status:** independent evaluator/design prepared and frozen before this reviewer accesses validation artifacts. No observed-result agreement or biological conclusion is claimed at this stage.

The final metadata/synthetic run passed checks for all2176 fixed classes,38338 distinct keys and20 architectures. It independently reconstructs categorical indicator ranks, complete codon-factor membership, nucleotide adjacency and connectivity. Least-squares projections agree with the protocol's conditional-mean projector on deterministic synthetic fields. Complete metric serialization/decomposition and max-min graph calculations were exercised on synthetic fields covering every architecture. The final Mathbox evidence record is `preparation_frozen/manifest.json`, status completed, and validates against the current source.

One useful metadata correction emerged: the eligible data represent162 development and154 validation central-AA strata, although the full possible-stratum assignment is200/199 plus quarantined VP. The failed reviewer assertion and source snapshot are retained. The upcoming influence denominators must use the represented eligible strata.

The scientific computations, exhaustive graph coverage and numerical comparison policy are specified in `design.md`. The independently authored `verifier.py` imports no author implementation. Its later evaluate phase requires root-supplied pinned validation artifacts. It checks all3 pseudocounts and computes every class's projection, residual/whole-field comparison, all45 normalization cross products plus the primary decomposition, influence contributions, and all9 scenario/raw/codon-additive barrier comparisons plus T/PT on all graphs.

The source protocol, inventory and training receipt are pinned in `freeze.json`. This reviewer inspected permitted training-export schemas and receipt, but did not open the workbook, RNA D, author class metrics or author validation summary. The root has indicated its own validation execution is complete; no numerical validation result was supplied to this reviewer. Evaluation and author comparison remain pending root notification after this frozen-ready report.

This is a separate review with inherited project context, not a blinded external review. All writes are confined to this owned review directory. No source/creator/repository/manuscript edits, commits, publications or biological interventions occurred. No independent biological calibration or statistical guarantee follows from this verifier's eventual numerical agreement.

## Phase 2: authorized independent evaluation and full comparison

**Verdict: computational agreement passes; the predefined scientific repeatability target fails.** Root supplied explicit permission and pinned validation artifacts after the verifier/design freeze. The frozen evaluator first computed and saved its complete independent output without reading author class metrics or the validation summary. Only then did the separately authored comparison adapter open those author outputs. No creator implementation was imported. The workbook was never opened.

The final comparison contains **3,352,385 checks and zero unresolved discrepancies**. It covers all2176 fixed classes at all3 pseudocounts, not a selected subset. Exact membership, key order, ranks, residual dimensions, strata, clique flags, training counts and export partitions agree. There are1124 development classes/20028vertices/162 represented strata and1052 validation classes/18310vertices/154 represented strata. No class was dropped for a null residual norm; the author and independent calculations both find no exactly zero primary training or validation class energies in this cohort.

The independent categorical indicator/lstsq projector agrees with the deposited training fields and the author's conditional-mean results. All293760 class/pseudocount normalization-pair decompositions, each with4 signed terms plus its direct total, agree. The primary averaged-training decomposition, both residual Gram arrays, RNA/DNA pairwise cosines and raw-count midrank correlations, pooled comparisons, losses, all948 leave-one-represented-stratum-out results, and influence fractions agree. All6 summaries preserve `target_failed`, and the earlier failed reference-protein pilot remains explicitly retained. The primary-summary alias agrees exactly with pseudocount0.5 validation.

Every selected graph was reconstructed independently from one-nucleotide Hamming adjacency. Max-min Floyd closure covered all endpoint pairs for all130560 field/graph calculations:9 scenario fields,9 codon projections, T and PT, for every class/pseudocount. All130560 author maxima and declared endpoint-pair counts agree within the frozen score-scale tolerance. Every62185 supplied nonempty maximum witness was checked against the independently reconstructed capacity matrix; the68375 empty witnesses correspond to zero reported maxima. Numerical-positive/tiny-positive flags, their declared machine-precision tolerance, held-out aliases and architecture-stratified summaries are consistent. This is exhaustive finite verification of the specified fields, not a proof about all possible biological landscapes.

| Quantity | Maximum absolute discrepancy |
| --- | ---: |
| Class residual cosine | 4.67e-15 |
| Class mean residual Gram entry | 1.78e-15 |
| Individual normalization signed term or total | 5.69e-14 |
| Pooled normalization signed term or total | 7.39e-13 |
| Stored T/PT/QT vector, Euclidean norm | 4.91e-14 |
| Graph maximum barrier | 3.11e-15 |
| Supplied barrier witness value | 3.78e-15 |
| Pooled primary residual cosine | 2.09e-17 |
| Leave-one-stratum cosine | 4.17e-17 |
| Pooled mean prediction losses | 0 |
| RNA/DNA raw-count rank correlations | 0 |

The complete category-level errors, locations, counts and empty discrepancy list are in `comparison_v3.json`. These differences are comfortably below the pre-outcome comparison budgets. No tolerance, endpoint, membership or frozen scientific formula was changed.

## Reproduced scientific result

The primary validation residual cosine is **0.07787615736918292**, far below the fixed0.5 target. Its leave-one-central-AA-stratum-out range is[0.06803624503075312,0.0919075805046242]. The fixed residual prediction has mean squared error0.11721404387352587 versus0.07934218820854008 for the zero-residual baseline, giving relative gain **-0.4773230549861418**. This means47.7% greater error under the specified descriptive loss. It does not estimate a biological noise ceiling.

| Pseudocount | Development residual cosine | Validation residual cosine | Validation gain against zero residual |
| --- | ---: | ---: | ---: |
| 0.25 | 0.05467078988646554 | 0.07923720520052298 | -0.48151748019005525 |
| 0.5 | 0.05867035361560697 | 0.07787615736918292 | -0.4773230549861418 |
| 1.0 | 0.06360366641102261 | 0.07763455283467703 | -0.4719660512387722 |

All six partition/pseudocount targets fail without rounding or reselection. The pseudocount checks and dependent stratum deletions do not rescue the primary result and are not independent confirmations or confidence intervals.

## Preserved adapter corrections and execution record

The frozen scientific evaluator remains byte-identical with SHA256 `4331bb0ed7349bfc896f59bdbdfb98b931818c6671b9fa908c8a261a5cdfa588`. Its independent output SHA256 is `b0a019cf1754598b672fba59af497ec9aefd1effabfa1515bb4cd15e695a2a82`. `evaluation_run_v1/manifest.json` records a completed, validated Mathbox run with unchanged inputs before/after. The evaluation took16.86seconds; the final comparison took16.33seconds on the recorded local runtime.

The first comparison adapter made two schema assumptions that failed: it compared raw Gram inner products to author class-mean products, producing117504 scale mismatches, and expected named library pairs where the author's metadata stores9 index pairs. The exact initial adapter is preserved as `compare_outputs_v1_snapshot.py`, with all117505 first-pass discrepancies in `comparison_v1.json`. The correction divides raw Gram entries and their propagated tolerance by the known class size and compares the identical indexed scenario order. It changes no mathematics, data or thresholds. `comparison_v2.json` then passes all2718531 original checks. The corresponding adapter is preserved as `compare_outputs_v2_snapshot.py`. The current adapter adds checks of all5 saved training log/residual library fields, raw training count arrays, target values/statuses, zero-energy counts, architecture counts, the primary-summary alias and retained-pilot flag; `comparison_v3.json` passes these additions as well. Nothing from an earlier failed run was deleted.

## Input identities and limits

| Artifact | SHA256 |
| --- | --- |
| Frozen protocol | 4ab4e85c14203c0ab8c983e8433c766f1566f373a335fd3fda4fb3981371c1aa |
| Decompressed inventory | 2e00b698dee109e357dc01435d5eacf8b341388f36463970454ade230664f350 |
| Training receipt | c6fb3309a93635efab38f1469516617c3f1b4f857c88b4fbfff91be21115bdef |
| Validation receipt | 5f7535bb58622cb656ee0cb7ed4109f8daef324ac470eb948ba40ffd9ab6757c |
| Author class metrics, compressed | f1c4c1fac09362edda68aa20b51c386784905837ca3eb9f89f256c7148e95491 |
| Validation count export, compressed | bcbf8989413b6fa5c9c1a3c945c75693f5b4a0e194a857754a41cc61848f09c9 |
| Author validation summary | cdbdd83d6861c54d1cd5094990cbeb306803349c63bdd0519caf391676c819f6 |

The receipt chain places training in commit3a10c6c74957bb7d3f559bf3067387926c1ed077 before author D access, with immutable train/validate copies supplied at3e88b7fdfdf4984edae2decd8541497182c4d69a. This audit independently checks the permitted exports and downstream computation, not raw FASTQ or workbook extraction. The original workbook identity remains inherited from the frozen protocol and receipts.

No computational correction is required in the author results reviewed here. The next action is faithful reporting of the failed target and negative fixed-prediction gain. Claims must remain conditional on these processed RNA/DNA library exports, the fixed pseudocounts and complete local classes in one truncated reporter. Local protein classes share most of the same construct and are not independent genes. Unresolved biological replicate identity, normalization contributions, latent nonlinear response and uncalibrated sampling error prevent causal interaction, selection, noise-ceiling or population-generalization claims. No new fit, threshold search, class selection, biological experiment, source edit, commit or publication was performed by this reviewer.
