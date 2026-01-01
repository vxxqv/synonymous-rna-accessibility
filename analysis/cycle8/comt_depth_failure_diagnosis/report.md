# COMT depth failure: source-contract diagnosis

**The unique post-collapse depth requirement was an added restriction, not a requirement of the source algorithm.** The frozen run correctly stopped under its declared contract, but that contract was an oversimplification of the source code. The failure is preserved. No pipeline repair, actual normalization, contrast, prediction, fit or model evaluation was performed in this diagnosis.

## Source evidence and precise behavior

The original notebook `comt_gate/analysis__integrative_MAVEs.Rmd`, SHA256 `ac267a0a5e1f07d3eadc30a2c3643a166349f0a67259ee31c649454229a869f0`, lines 2148-2150, calls `postprocess_vcf_summary_dt` with sample/ROI/input/fraction/quadrant/variant keys. The held `comt_pilot/summarization_utils.r`, SHA256 `f0a741b087e9e06e6d32cd947f2faf2dcedd59272d3e0a87921de118bad740bd`, lines 196-241, sorts by those keys, takes the median of numeric per-base fields, takes the first nonnumeric field, and merges the results by the keys. DP is among the numeric fields.

The notebook then uses:

```r
# lines 2177-2178
.(CAO_sum=sum(CAO), DP=unique(DP)),
    by=.(ROI, Input, Fraction, Quadrant, Sample, POS)
# line 2180
WT_CAO := DP - CAO_sum
# lines 2183-2187
WT_key := paste(Sample, POS, sep=":")
match(comt_all_filt_postprocess_dt[,WT_key], WT_key)
```

`unique(DP)` may return multiple values. In the reviewed installed data.table runtime, that group expression emits multiple denominator rows and repeats the scalar CAO sum. `match` selects the first matching `Sample:POS` row for every allele at that start position. There is no assertion that only one DP exists. The source does not select each allele's own DP, average depths, or explicitly select their maximum. The first DP is determined by the source-collapsed table ordering. The helper sorts by raw variant ID within the fixed sample metadata; the first record is therefore not a random input-row choice in the tested fixed runtime.

NC subtraction occurs later: notebook lines 2194-2211 compute per-raw-variant control medians, join them to other samples, use zero for absent controls, subtract normalized frequencies and retain positive corrected frequencies. This later operation cannot explain or repair the current pre-subtraction uniqueness failure.

## Complete observed depth audit

The diagnostic verifies the frozen extraction/status seals and reads only the permitted identity columns, DP and Boolean frequency/rejection flags. CAO, CAF and classifier covariates were not numerically decoded; actual ratios, WT values and contrasts were not calculated. All existing status masks were used unchanged.

| Diagnostic | Count |
| --- | ---: |
| Extracted rows inspected for identity and DP | 327,584 |
| Rows retained by the existing frozen flags | 252,872 |
| Collapsed raw alleles | 111,448 |
| Retained sample/start-position groups | 5,535 |
| Groups with more than one median DP | 1,091 |
| Samples containing a conflict | 61 |
| Start coordinates containing a conflict | 155 |
| Conflict groups containing a selected target | 364 |
| Direct target/sample denominator records in conflicts | 740 |
| Distinct selected target IDs affected | 116 of 125 |
| Affected target IDs by type | 105 SNV; 11 non-SNV |
| Matching selected-target NC records in conflicts | 0 |

Conflict counts by assay role are DNA 176, RNA 348, low-bin FACS 132, heavy 206, mid 219 and NC 10. ROI1 has 425 conflicts and ROI2 has 666. These are data-contract counts, not biological effect or prediction results.

Every one of the 1,091 conflicts contains SNV and multinucleotide alleles. Within a retained raw allele its per-base DP records are constant, so median aggregation is **not** the observed cause of the discrepancy. SNVs at a common start agree on DP; at least one multinucleotide allele has a lower DP, and no multinucleotide allele exceeds the SNV depth in these groups. Thus the deposited DP is associated with the allele's counting object, rather than universally equal for all alleles sharing a start coordinate. The upstream caller's exact reason for the lower multinucleotide depth is not established by the inspected helper and notebook; a read-span/eligibility explanation remains a hypothesis, not a verified cause.

For all 1,091 conflicts, the first raw variant ID under the source key ordering is an SNV, so the literal first-match source rule would choose the shared SNV depth here. This empirical equivalence does not justify replacing the source rule by a universal maximum-depth rule for other inputs. No actual WT denominator or normalized outcome was chosen in this audit.

`all_conflicts.csv` includes every conflict. `all_conflict_variant_lineage.json.gz` records every participating allele's raw record IDs, positions and depths, collapsed median, target identity and first sorted variant. `affected_target_sample_rows.csv` gives all 740 affected target/sample records. `receipt.json` lists every affected target ID, input hashes and the original failure hashes.

## Minimal invented demonstration

`source_rule_invented.R` uses the unchanged author helper and the literal notebook WT expressions in the pinned task-local runtime. Two SNPs have DP100 and CAO10/20; one two-base allele has DP90 and CAO5 on both per-base records. Collapse produces three alleles, so CAO sums to35. The source emits DP100/90 and WT65/55, then assigns the first value65 to every allele. The Python uniqueness guard would reject the same object. Reversing raw input records before the source helper preserves the relevant collapsed IDs/DP/CAO. Deliberately reversing the already-collapsed table changes the first-match WT to55, demonstrating the ordering dependence of the source expression itself.

The passing invented output is in `invented_R_v2/`. An initial test used `identical` on the entire R data.table objects and failed before the WT assertions. Its source and failure are retained in `history_v1/`. The corrected test compares the relevant collapsed VAR_ID/DP/CAO fields using `all.equal`; no source helper or biological data was changed. The distinction between finite demonstration and source semantics remains explicit.

## Successor options

1. **Literal source-rule successor in the already fixed new runtime.** An explicit amendment can replace the added uniqueness restriction by source-equivalent ordered `unique` plus first-match selection, while retaining full denominator lineage and every original filtering/target rule. This is the narrowest source-fidelity correction. It needs independent comparison against the R source on varied invented orderings, tied depths, multiple distinct depths and both variant types before a new freeze. It remains a new-runtime reconstruction, not recovery of the author's historical fitted state.
2. **A different scientifically motivated depth definition.** Choosing a canonical per-base depth, an allele-specific depth, a median, minimum, maximum or other coverage rule would define a different normalization. It needs an upstream count/depth contract establishing the relevant read universe and allele partition. It cannot be adopted merely because it eliminates the failure, because one choice matches the observed first value, or because it improves eventual prediction.
3. **Close this raw-count branch.** If source ordering cannot be authenticated sufficiently or the source-defined denominator is unacceptable for the intended proxy question, preserve the failed reconstruction and retain only the data-contract limitation. Do not silently drop the 740 records, affected alleles, samples or sites.

Option1 is **not merely a file-format amendment**: it changes the admissible numerical domain and the executable denominator-selection rule relative to the frozen contract. It is a source-fidelity normalization-contract amendment that must record this observed diagnostic exposure, preserve the failed run, receive fresh review and be frozen before further normalization. Option2 is a more substantive measurement-definition amendment. Neither option establishes that source `WT_CAO` is an independently measured physical wild-type count; DP minus a sum of variant counts still depends on the caller's read-universe and exclusivity assumptions.

## Preservation and exposure boundary

All original extraction, RF, protocol, kernel and implementation files remain unchanged. The first normalization failure is copied byte-for-byte under `preserved_failure/` and retains its original location. No source RF was rerun. No source-status flag, threshold, subset, fold, target or baseline changed. The observed information is limited to identities, the frozen inclusion mask, allele DP and its conflict structure. There are no observed q, NC-adjusted phenotype, R/T/L, prediction or model-performance results in this package.
