# COMT depth-only failure diagnosis

Verdict: the frozen unique-depth stop is independently reproduced. The original source executable uses a different rule. Restoring that rule is a measurement-definition amendment after the failed normalization, not a file-format repair. No actual allele count, frequency, reconstructed WT, normalized ratio, log contrast, prediction or performance result was decoded by this audit.

## Independently reconstructed membership and depths

The separate depth_only_probe.py was hash-sealed before actual DP-field decoding. Its binary selective reader decodes only record/sample/reference/variant/base identities and DP; only the four named RF row-status fields are decoded from the mask files. A poisoned invented CAO/CAF fixture and -O/-OO refusals passed first. The source TSV files were byte-hashed against the two authorized stage seals before selected-field parsing. Hashing their bytes is not decoding their excluded columns. No production scientific function was read or imported.

| Quantity | Independent result |
|---|---:|
| Extracted records with matched RF masks | 327,584 |
| Retained per-base records | 252,872 |
| Retained collapsed raw alleles | 111,448 |
| Collapsed alleles with variable per-base DP | 0 |
| Retained sample/start groups | 5,535 |
| Groups with more than one collapsed allele depth | 1,091 |
| Distinct samples / starts in conflicting groups | 61 / 155 |
| Conflicting groups containing a target allele | 364 |
| Direct target/sample records in conflicting groups | 740 |
| Distinct impacted target IDs | 116 of 125 |
| Impacted SNV / non-SNV target IDs | 105 / 11 |

The 740 records and 364 groups are unchanged if negative-control samples are excluded. They count direct records whose own position's denominator is implicated. They are not completed valid measurements, common model coverage, 116 independently affected experiments, or an estimate of the eventual performance impact. A nine-ID complement is not a new eligible cohort.

All 1,091 conflicts have a single common raw-SNV depth. Every MNP depth is at or below it and at least one MNP depth is strictly below it. In 322 groups some additional MNPs have equal depth; in 769 every MNP depth is strictly below it. No MNP is above it. Thus the statement "all MNPs are strictly lower in all conflicts" would be false; the source agent's actual weaker wording, some lower and none higher, is confirmed. The median collapse itself is not causing these differences: no retained raw allele has variable per-base DP.

The source-lexical first raw VAR_ID within every conflicting sample/start group is a one-base SNV at the common SNV depth. This is relevant to the declared C-collation reconstruction and source helper ordering. It does not authenticate an unspecified historical runtime/classifier state.

## Source rule and why it matters

The pinned original notebook work/research_cycle8/comt_gate/analysis__integrative_MAVEs.Rmd, lines 2148-2150, invokes postprocessing with Sample, ROI, Input, Fraction, Quadrant and VAR_ID keys. In the held original summarization_utils.r, line 144 creates the raw VAR_ID factor. Lines 196-241 key/sort, collapse numeric columns by median and merge the keyed tables. The extra assay labels are fixed within each selected sample in the source metadata contract, so the source grouping distinction does not create an alternative DP value within a sample/start group.

The notebook then executes:

- Lines 2177-2178: CAO_sum=sum(CAO), DP=unique(DP) within the assay/sample/start group. Multiple distinct DP values are permitted by that expression.
- Line 2180: each resulting candidate has WT_CAO=DP-CAO_sum.
- Lines 2183-2187: build Sample:POS keys and use match, which selects the first duplicate denominator row.
- Line 2189: compute (CAO+0.5)/(WT_CAO+0.5).

This differs from the frozen protocol's requirement that collapsed depths agree. A new implementation cannot call the frozen failure an input-format issue or silently bypass it. The declared pipeline correctly stopped under its specified stricter assumption, but that assumption was not an exact description of the executable source normalization on these data.

The base-R invented test source_match_synthetic.R uses depths 100,100,90 and allele counts 10,20,5. unique(DP) yields 100,90, the candidate WT counts are 65,55, and match assigns 65 to every allele. Reversing the candidate denominator table assigns 55. With source-ordered SNV depths 80,80 and MNP depth 100, first-match instead assigns WT45, whereas a maximum-depth rule would assign WT65. Therefore "take the maximum" is not a generally equivalent replacement for first-match, even though the two depth choices coincide at all 1,091 present conflicts. No actual WT or q was evaluated to establish this structural fact.

All assertions in this separate base-R synthetic test passed. Its startup emitted four C.UTF-8 locale-setting warnings, retained in the review history; the script explicitly set LC_COLLATE to C before its ASCII ordering assertion. This is not represented as a warning-free source replay or a rerun of the complete author helper.

## Consequence and stopping rule

The smallest defensible successor is an explicitly versioned, reviewed source-rule amendment preserving the original failed normalization, source sort/key semantics, first-match denominator choice, and the still-required negative-WT/domain checks. It would be a post-extraction, post-failure measurement amendment before any contrasts or predictive results, not an unchanged preregistered primary or confirmed historical reproduction. Reuse of the old source-state masks must remain explicit. A new broad maximum, minimum, per-allele denominator, allele-class exclusion or favorable target subset would be an additional method, not this source rule.

The present exposure cannot establish that revised WT counts are nonnegative, that negative-control subtraction yields available logs, that model coverage reaches its threshold, that either FACS map is correct, or that heavy load improves prediction. No physical explanation for the different read depths follows. In particular, the depth pattern is not proof of a biological MNP effect or a valid shared physical read denominator.

This closes only the depth-assumption failure. A future source-rule replay must stop on a mismatch or any invalid reconstructed denominator, retain all unavailable records and existing failed hurdles, and receive fresh explicit authorization. This audit performs no such replay and changes no production source, frozen reference, threshold or model.

## Reproduction

From the workspace, with an available Python interpreter and the exact hash-pinned source artifacts:

```text
python depth_only_probe.py --self-test
python depth_only_probe.py --extract work/research_cycle8/comt_results/actual_v1/extract --rf work/research_cycle8/comt_results/actual_v1/rf --contract work/research_cycle8/comt_pilot/implementation_v2/selective_contract.json --out NEW_DP_ONLY_DIRECTORY
python check_mnp_relation.py --input NEW_DP_ONLY_DIRECTORY/conflicting_groups.csv --out NEW_MNP_RELATION.json
Rscript --vanilla source_match_synthetic.R NEW_SYNTHETIC_TEXT_RECEIPT
```

Script paths are relative to this failure_audit directory; input paths in the example are relative to the workspace. All destinations must be new. The probe does not read the extraction's simulated training files, RF prediction probabilities/model objects, normalization outputs or any downstream result. Its pre-access seal, actual-run receipts, complete conflict and direct-target tables, and later relation check remain separate and unchanged.
