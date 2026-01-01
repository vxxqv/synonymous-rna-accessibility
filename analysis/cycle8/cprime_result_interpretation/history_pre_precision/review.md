# cPRIME B1 scientific interpretation review

**Verdict: a defensible finite local-assay sensitivity result, with sensitivity-fragile ordinal exclusions. No isolated intended-edit effect, edit-marker interaction, calibrated population statement or new processing mechanism is identified.** This review addresses interpretation; the separately assigned independent numerical replay is not replaced by the checks here.

## Complete result and widening sensitivity

The saved result retains all 30 fixed intended labels and 66 local marker states, both cells, 132 state/cell records and 396 source rows. All recorded state intervals are complete; no missing state was dropped. The 60 group/cell entries are dependent contexts from one study, not 60 independent experiments or a random prevalence sample.

| Fixed envelope factor gamma | Lower spread >1 | Lower spread >=2 | Upper spread <=2 | Exons excluding one common order |
| --- | ---: | ---: | ---: | --- |
| 1 | 19 | 7 | 13 | 17,27 |
| 1.25 | 7 | 5 | 4 | none |
| 1.5 | 6 | 4 | 0 | none |
| 2 | 4 | 2 | 0 | none |

These are the protocol's reported exact finite-set quantities, not p-values. At gamma=1, 19 contexts cannot share one marker-independent response while each response remains inside its three-row empirical range; 41 permit equality in some completion. Seven have an unavoidable factor-at-least-two spread under that assumption. At gamma=2 four still exclude equality and two retain a factor-at-least-two lower bound. That sensitivity survives the declared enlargement, not every measurement error or biological preparation. The factor two comparison is an operational scale choice.

The 13 primary upper-bound certificates say every completion in the primary box has spread at most two; they are not demonstrations of marker neutrality. At gamma>=1.5 there are no such certificates. This does **not** say every context truly exceeds twofold: a large upper bound often means uncertainty, not an unavoidable difference. Lower- and upper-bound counts must not be treated as a binary effect/no-effect classifier. The entire fixed 1.25,1.5,2,4 comparison grid remains in evidence.json and the primary result.

The empirical box permits independently chosen responses inside their preparation envelopes. Its completions need not correspond to one actual jointly observed preparation. The larger boxes are prescribed possibilities, not calibrated confidence sets, measured correction factors or estimates of how likely an artifact is. No optimal radius or new threshold was sought in this review.

## Ordinal certificates and their actual scope

All order comparisons are confined to one of exons 17,25,27,39,42. No cross-exon order or numerical scale is tested. Each recorded cycle edge matches a strict upper-versus-lower inequality in its named cell's state interval. Ties are not converted into edges.

The exon-17 two-cycle compares two different marker states of the **same** intended label chr17-31225238-G-T. KBM7 forces one order and K562 the reverse. The exon-27 reported three-cycle spans **three different intended labels**, not three markers of one intended edit. It therefore excludes a common ordering of the exon-level local states, not specifically a within-edit marker reversal. Both certificates contradict the stated unrestricted scalar score plus separate nondecreasing cell maps, conditional on the primary empirical intervals and source identity.

Both exons become feasible already at gamma=1.25; all five are feasible at every wider factor. Say that these primary exclusions do not persist under the first declared enlargement. Do not say that 25% measurement noise was measured, that the correct common score was recovered, that the model became probable, or that an additive model now fits. Widening both ends allows a pairwise separation ratio to change by gamma squared. Feasibility is merely an existential monotone completion of that enlarged box.

## The contrasts do not require interactions

Within an intended label, different marker main effects alone can produce every observed kind of spread: Y_em=mu+alpha_e+beta_m with eta_em=0 already gives marker dependence. Without isolated edit and matched control states, alpha and eta are not separately determined. The excluded arrayed panel does not repair that deficiency.

Even the same-intended-label cell reversal needs no edit-marker interaction. In a synthetic two-marker example, marker states (x,y)=(1,0) and (0,1) have log proxy x+2y in cell A and 2x+y in cell B. Exponentiating gives positive responses 2 versus 4 in A and 4 versus 2 in B. Both cells have strictly additive, cell-specific unary marker effects and zero edit-marker interaction, yet no common scalar marker order exists. This is an analytic counterexample to an interaction interpretation, not a fitted explanation of the data. The exon-27 cross-label cycle allows additional context-specific intended main effects as well.

There is no measured SNV path or complete mutant graph in this pilot. Switching marker states can change different nucleotide positions. Neither a barrier nor preservation of a processed protein or full haplotype follows from these tables.

## Most damaging alternative explanation

One common positive scaling of every state in an exon/cell cancels from marker fold spreads and preserves all order constraints. A generic claim that such a single normalization constant explains away the result is mathematically wrong.

Preparation- or state-dependent normalization and detectability are different. With the source count convention, a row ratio has the form

    q = [(RNA_count+1)/(DNA_count+1)] * [(DNA_total+1)/(RNA_total+1)].

The totals factor can differ across anonymous preparations; the processed table omits the row-to-library/replicate crosswalk and source count totals needed to retain that structure. Preparation-specific scaling may change envelope overlap even though a single fixed exon/cell factor cannot. Unknown low-count pseudocount influence cannot be bounded from processed RPM and score agreement alone. Source reproduction checks arithmetic of the deposited estimator, not its calibration to absolute biological abundance.

More seriously, the local search word may represent multiple unphased realised backgrounds, while RNA recovery requires that word and primer sites in a nearcanonical selected product. State/cell-specific changes in isoform visibility, amplification, sampling or background mixture can change q without changing an isolated intended-edit effect or total RNA. In an idealized positive measurement model q=a*rho*t, even constant biological t is compatible with any positive proxy vector if state/context-specific rho is unrestricted. This is a nonidentifiability countermodel, not evidence that such bias actually occurred or an exact raw-count likelihood for the pseudocounted data. The strongest surviving statement remains about local-sequence-conditioned assay representation.

Shared preparation, denominator and recovery effects make the 60 entries dependent. An empirical range over three source rows is not a confidence interval; their lack of explicit pairing does not supply independence. Conservative statistical control could use authenticated marginal coverage without guessed paired covariance, but no such coverage or systematic-recovery bound was authenticated here.

## Closest prior and bounded increment

The original cPRIME marker section already demonstrates that combined edits and markers can differ from edit-only responses, investigates marker alternatives, discusses splice-motif explanations and warns against interpreting marker sequences as neutral. The selected arrayed examples are exposed prior, not independent validation of this audit. The B0 T-marker motif discrepancy remains an unresolved local identity record and was excluded from this pilot; the current result neither repairs it nor demonstrates that the original physical experiment was wrong. [El Kassem et al., cPRIME, Nature Communications (2026), DOI10.1038/s41467-026-78344-9](https://doi.org/10.1038/s41467-026-78344-9).

The incremental result is the complete prespecified 30-label, all-marker finite-set audit, with sharp lower/upper dispersion, explicit widening sensitivity and exact common-order exclusion certificates. It is not discovery of marker non-neutrality, a new generic interval theorem, a processing mechanism or confirmation of original predicted-splicing examples. Different scores and models cannot substitute for the missing observation contract.

## Failed or inconclusive stronger conclusions

- Universal finite-marker equality is excluded only for the stated subset and envelope; the remaining equality-compatible contexts are non-exclusions, not neutral markers.
- All primary ordinal exclusions fail to persist at the first frozen widening. No shared additive representation is certified by the feasible wider cases.
- No biological significance threshold, calibrated confidence, population frequency, clinical effect, isolated intended-edit sign, interaction magnitude, total-RNA effect or protein/isoform preservation is identified.
- Cross-cell results remain within-study context comparisons. No new-study or independent-culture replication is established by this analysis.

## Single smallest missing source record

The closest missing source object is an authoritative **396-row-to-library manifest**: processed row identity mapped to its exact paired RNA/DNA library accessions, cell/exon, biological/preparation/technical replicate identities and the RPM denominator totals. That one mapping would expose which offsets and denominators are truly shared and make the dependence structure auditable without inventing pairing. Its absence does not invalidate the descriptive arithmetic. Recovering it would still not itself calibrate the widening factors, remove state-specific recovery, identify isoforms or isolate an intended edit. Those stronger conclusions remain separate gaps rather than reasons to fabricate a successful pilot.

## Provenance and limits

Authorized result SHA256 eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf; freeze SHA256 58ef31aeed9a23aeb13a8ab54f3c74d4224b410aad71c756a7d3a612623b4363. The root reports public freeze commit b6f83e550a8e6f7bbf7901b9e9d69e01956eaac5. This review checked current frozen file hashes and saved result linkage, not remote historical access order. It did not rerun the scientific analysis or independent numerical replay. Sources were the frozen protocol, saved result/receipt, eligibility metadata, B0 source contract and cached original primary marker prose. No source workbook, new raw counts, COMT outcomes, model, fitted nuisance factor, new threshold, favorable subset or external contact entered this review.
