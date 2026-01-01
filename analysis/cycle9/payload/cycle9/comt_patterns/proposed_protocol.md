# COMT local readout-order investigation: proposed bounded protocol

**Status: ready for root review, no new observed comparison computed.** This is an exploratory follow-up to exposed Cycle 8 data. It cannot become prospective validation by making a new split. The frozen prediction failure, source-depth failure and every subsequent amendment remain unchanged. The question is a local ordinal limitation of these measured proxies, not a rescue of assay substitution.

## Question and possible additional contribution

Among all authenticated selected synonymous alternatives at the same COMT codon, is there an allele pair for which one allele has higher RNA representation and lower oriented low-fluorescence-bin contrast under every represented preparation combination and both source-supported FACS grouping maps?

The additional empirical object would be a **complete inventory of preparation-uniform, same-codon allele-order reversals**. A strict reversal excludes a common nondecreasing RNA-to-proxy response for that pair in each represented scenario. Hoskins et al. already reported discordant expression layers and site-specific effects. Cycle 8 already reported four variants at three sites with opposite variant-versus-WT proxy signs. Neither finding is to be relabeled as the new result. The comparison here is alternative versus alternative, retaining both alleles and their shared assay preparations.

The mathematical implication from reversed order to failure of a nondecreasing response is elementary. This is not a new theorem. If no new local ordering survives, report that informative limitation and stop. Do not move to global cross-site pairs, change the sign threshold, select a favorable map or add a predictor.

## Fixed identity and coverage universe

Use the byte-pinned Cycle 8 v3 source contract and its 125 selected synonymous nucleotide identifiers at 53 sites. All were translated against the authored fusion in the previous source audit. This does not authenticate all endogenous isoforms, full clone haplotypes, or a new independent biological context.

Metadata alone gives 134 unordered alternative-allele pairs across 32 sites with at least two selected alleles. Keep the 21 singleton sites in the coverage report as not informative for an alternative-versus-alternative comparison. There are 98 pairs with edit-string distance one and 36 with distance two; these are descriptive labels, not evidence of evolutionary routes. The fixed 113-SNV/52-site subset contains 97 pairs across 31 sites and 21 singleton sites. It is reported from the same calculations without refitting or replacing the full universe. Preserve ROI and every site's pair count.

Order the two alleles lexically by UTF-8 nucleotide identifier solely for an immutable sign convention. Never choose orientation from outcomes. `metadata_inventory.json` contains every pair and metadata denominator. Missingness may make a pair unresolved but never remove it from the denominator.

## Existing observations and exact estimands

Read only the sealed `actual_v3/normalize/contrasts.json` after approval. Its assay-group medians and technical sample identities were already used and independently verified. No count extraction, RF fitting, negative-control outcome retrieval, ComBat, ALDEx2 or normalization rerun is needed. Do not read any new source-control value.

For allele v, map m, assay A and declared group g, let x(A,m,g,v) be the existing median of finite technical natural-log source-normalized representations. These are not absolute molecule amounts. Let the assays D, R and F denote gDNA, total RNA and the source low-fluorescence P3 population, respectively. Four groups per assay are declared. Do not call them paired cross-assay replicates.

For the ordered pair (a,b), define the same-assay paired difference

`z(A,m,g;a,b) = x(A,m,g,a) - x(A,m,g,b)`.

Preserve the source's median for each allele, then subtract. The median of technical pair differences is a different estimand and is not substituted. Both members must have the same declared technical sample list and the same nonmissing contribution mask in every needed assay/group. If a mask differs, retain the pair as unresolved with the exact mismatch; do not select their common subset or reaggregate it.

For every i,j,k in {1,2,3,4}, define

`DeltaR(m,i,j) = z(R,m,i) - z(D,m,j)`

`DeltaT(m,j,k) = z(D,m,j) - z(F,m,k)`.

Thus each map has 64 declared joint choices, with the same gDNA group j used in both coordinates. RNA and FACS preparations are not matched by numeric labels. Every source-supported cross-assay combination is retained. The maps are alternative interpretations of shared data, not independent experiments.

Positive DeltaR means the a-versus-b contrast is greater for total RNA relative to gDNA. Positive DeltaT means a is less enriched in the low-bin proxy after its gDNA contrast is accounted for. It does not mean a known increase in mean protein abundance. Neither polysome load nor the target is translation flux or fitness. Heavy and mid-load fields remain in the original frozen result; no new load model or ranking is introduced here.

## Primary exact classification

Convert each finite stored binary float to its exact rational value for arithmetic. This is exact accounting of stored numerical inputs, not exact biological logs. Unknowns remain unknown; NaN or invalid schema stops. No extra numerical or biological tolerance is introduced. An equality or zero is not a strict direction.

For each map save the four same-assay differences, all 64 joint coordinates, the minimum/maximum DeltaR and DeltaT, and their attaining source group identities. Retain tied attaining identities. Combined bounds take the extrema over both maps. These extrema also describe the Cartesian hull of the same-assay difference ranges, with shared D retained in joint records.

The main positive predicate is either:

1. minimum DeltaR over both maps is strictly positive AND maximum DeltaT is strictly negative; or
2. maximum DeltaR is strictly negative AND minimum DeltaT is strictly positive.

The same orientation must hold under both maps. The quantifier is one fixed pair and one fixed direction, for every represented map/preparation combination. It is not existence of a different favorable pair in each scenario. This can be checked both by the extrema formula and by an independent enumeration of all 128 joint records.

Also report same-direction strict separation in both coordinates, per-map classifications, exact ties, missing cases and otherwise unresolved comparisons. These categories form an explicit exhaustive table. Call the construction a strict-separation partial order, not a complete Pareto frontier: weak dominance with ties and all possible biological objectives are not claimed. Counts of reversals among recorded centers may be reported only as descriptive context, with the uniform result kept separate.

Do not infer existence of a common monotone response from failure to find a uniform reversal. That failure does not settle weaker scenario-dependent reversals, collective ordering constraints, or an unrestricted joint model. Do not fit those further models in this pilot.

## Sensitivity and dependence

For a uniform reversal let m be the smaller distance of its two combined contrast ranges from zero. Under independent additive perturbations of magnitude at most delta to each allele's assay-group log value, a pair contrast changes by at most 4delta. Its strict reversal is preserved exactly for delta < m/4; a tie is permitted at equality. Save this deterministic margin for every pair, with zero for a noncertified pair and null for missing input. Do not choose a delta cutoff to declare a biological effect. This Cartesian allowance is not a confidence interval, measurement-error model or bound on systematic source normalization bias.

Report the complete 134-pair table, fixed SNV flags, all 53 sites including the 21 noninformative singleton sites, per-site fractions and pair concentration. Provide both pooled pair counts and equal-site descriptive summaries over the 32 informative sites, with explicit denominators. Neither pairs nor sites are independent experimental replicates; alleles, preparations, library construction, WT denominators and negative controls are shared. No confidence interval, p value or independent-study claim is warranted.

A uniform result cannot be created by selecting a single favorable preparation from this represented set. List all limiting preparations and all nonuniform cases. Do not remove a limiting preparation and promote the resulting classification as the primary result.

## Mandatory normalization challenge

The defining algebra is `DeltaR = zR - zD`, `DeltaT = zD - zF`. A synthetic field with zR=zF=0 and zD=-1 gives DeltaR=1 and DeltaT=-1 in every preparation. It can therefore satisfy the main predicate without a difference in the two assay-specific numerator log representations. A shared DNA contrast alone can create the opposing proxy signs. Even a large uniform reversal is not an intrinsic RNA/protein trade-off, kinetic compensation or an evolutionary fitness trade-off.

Within a source sample, equal position-WT denominator Z cancels from an uncorrected log ratio between two alleles only before variant-specific negative-control subtraction: qv=av/Z. With source corrections cv and cw, the log ratio becomes `log(av-Z*cv)-log(aw-Z*cw)`, so dependence on Z generally remains. Same codon, same edit start or common library does not prove normalization independence. Source first-depth selection, position-specific denominators, variant-specific control subtraction, pseudocounts and RF ascertainment remain conditioned assumptions. Do not test a new denominator or omit negative controls after seeing a pattern.

Save the existing descriptive numerator differences zR,zD,zF beside every pair so readers can see the algebra. This is a decomposition, not a counterfactual re-normalization. Identifying a physical trade-off requires a newly authenticated, calibrated and appropriately replicated assay or control. That stronger branch is outside the current authorization.

## Synthetic checks before observed computation

- Uniform opposing ranks across both maps: correct orientation and all 128 joint cases.
- Uniform same-direction ranks: no reversal.
- One differing preparation or map: not uniform; retain the offending identity.
- Zero/touching ranges: not strict; no tolerance rescue.
- Missing group or unequal technical contribution mask: unresolved, not dropped.
- Swapping a and b: reverse signs, retain classification and margin.
- Exact common within-assay offsets for both alleles: pair differences invariant.
- The shared-D-only construction above: predicate positive, mechanistic interpretation explicitly invalid.
- An invented unequal-NC example: common Z does not cancel after source subtraction.
- Perturbation boundary: strict classification below m/4 and permitted tie at equality, independently enumerated.

The independent verifier should reconstruct pair classification directly from the finite group records without importing the primary classification function. Retain complete failures. Implementation and synthetic receipt must be frozen before opening the observed comparison outputs.

## Decision and stopping rules

The smallest adequate run is the complete 134-pair comparison with both maps and all declared group combinations. No new fit, sequence covariate, source-control measurement, threshold search or held-out split is needed. If at least one uniform local reversal exists, it supports a named-pair ordinal limitation of these source-normalized proxies under the declared finite uncertainty set. It is exploratory same-source evidence, not an independently validated biological mechanism. Compare every surviving example with the original paper before claiming an additional empirical observation.

If only WT comparisons already reported, nonuniform cases, or prior-published equivalent local comparisons remain, record no distinct new biological result. If no pair is uniform, report the exact coverage and margins without replacing this test. If identities, source hashes, technical sample alignment or source implementation differ from the pinned contract, stop and retain the discrepancy. Any new source control, normalization definition or follow-up model requires a separate question and approval.

No result-dependent expansion is authorized by this proposal. Root review is the next step.
