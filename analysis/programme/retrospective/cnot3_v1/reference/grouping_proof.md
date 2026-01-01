# Proposed near-identical-protein grouping, not yet applied

This is a proposed deterministic leakage precaution for the eventual fixed eligible gene cohort. It is not an executed partition, a biological family annotation or a claim that genes become independent. Root must freeze the cohort, normalization, threshold, resource gate and split assignment before outcomes or grouping results guide choices.

For two nonempty complete representative proteins, define `d(s,t)` as full-length Levenshtein distance with unit substitution, insertion and deletion costs, and define normalized distance as `d/max(len(s),len(t))`. The proposed edge rule is exactly `10*d <= max(len(s),len(t))`. The maximum-length denominator is an explicit proposal, not an implied convention. Combine exact identities and these edges, and keep every connected component in one partition. Connected components deliberately include chains whose distant endpoints can exceed10% distance. No pair meeting the declared edge rule may cross partitions. Near-identical proteins outside the eligible cohort need not be processed, but an optional broader family grouping would require separately frozen evidence.

## Lossless candidate bound

Use the multiset of all overlapping5-mers of the protein. Let its cardinalities be `N_s=max(0,|s|-4)` and `N_t=max(0,|t|-4)` and let `I=sum_w min(count_s(w),count_t(w))`. Put `D=floor(max(|s|,|t|)/10)`. A pair can pass only if:

```text
abs(|s|-|t|) <= D
max(N_s,N_t) - I <= 5*D.
```

Proof: one substitution deletes and introduces at most5 overlapping5-mers. One insertion deletes at most4 and introduces at most5; deletion reverses this. Boundaries only decrease these numbers. For multisets `A,B`, the quantity `max(|A|-|A intersect B|, |B|-|A intersect B|)` is the maximum of the two sums of positive multiplicity differences. Each directional sum obeys the triangle inequality, so their maximum does too. Along an edit sequence of length `d`, the displayed multiset deficit is therefore at most `5*d`. This proves a necessary bound, not sufficiency. Using unique sets instead of multisets would invalidate the guarantee, especially for repeated amino acids.

An inverted index from5-mer to unique protein records can accumulate `I` by adding the minimum multiplicity for each shared5-mer. First collapse exact protein identities while retaining every member gene. For lengths below10, the distance budget is zero and identities are the only eligible edges. For maximum length at least10, `L-4-5*floor(L/10)>0`, so any passing pair must share at least one5-mer and is present in the inverted-index candidate union. Thus this candidate generation has no false negatives for the declared edge rule. Apply the length and multiset bounds, then evaluate the exact distance for every surviving candidate. Deterministic sorted IDs and sorted5-mers remove dependence on traversal order.

## Exact implementation option and resources

`grouping_proposal.py` demonstrates the exact scoring primitive using the installed C-backed Biopython `PairwiseAligner`: global mode, match0, mismatch-1, all gap-open and gap-extension scores-1. A run of `q` gap characters costs exactly `q`. Negating its optimal global score is the unit-cost edit distance. Only the score is requested. All penalties are fixed; there is no local alignment, free terminal gap or protein substitution matrix. Integer scores are required. `test_grouping_proposal.py` independently compares this score with a separately implemented integer dynamic program.

Before any alignment, a frozen resource gate can bound index work by `sum_w choose(posting_count(w),2)` and subsequently bound exact dynamic-program work by `sum_(s,t) |s|*|t|` over retained unique-protein candidates. Candidate bookkeeping can still be quadratic for common/repetitive5-mers. The bound does not promise subquadratic worst-case runtime. Each exact global comparison has an `O(|s||t|)` time upper bound; no stronger implementation-specific memory guarantee is needed here. Root should choose explicit count/work/wall-time caps before executing on the eventual cohort. Exceeding a cap must stop the analysis and report an unmet resource gate. It must not drop candidate edges, raise a similarity threshold, split an oversized component or select a smaller favorable cohort. No numeric cap has been silently selected in this proposal.

## Checks and interpretation limits

Synthetic checks passed for3,969 independent dynamic-program/C-backed comparisons,11,917 multiset bound checks and10 threshold/boundary fixtures, including substitutions at exactly10%, insertions/deletions, repeated amino acids and short identity-only cases. No real reference protein pair was aligned and no grouping or partition was applied. The finite tests supplement the lossless bound proof; they do not prove every library implementation detail or establish biological relatedness.

The existing exact-protein identity option groups173 nontrivial components involving408 of the19,567 selected genes. Its weaker guarantee remains explicitly labeled. A10% full-length edit rule still misses more distant paralogs, domain-level homology and regulatory-family dependence. Passing a group-aware holdout would support only a declared predictive generalization on this reporter/reference/context, not causal codon effects, independent biological units or broad family-independent prediction.
