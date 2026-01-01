# Gill multiclass sequence eligibility and reproducibility proposal

**Decision: a substantial finite, within-reporter reproducibility analysis is eligible for a new freeze. Independent biological replication, a noise ceiling and transfer across genes remain unidentifiable. No nonreference outcomes have been opened.**

## Exact denominator

The authenticated forward orientation, pMG23 reporter and final amplicon identify the same eight-base variable interval at CDS offset 467, affecting codons 156-159 of a 313-residue truncated aprE reporter. This audit independently reparsed the construct, compared a literal codon table with Biopython for all 65,536 strings, and compared all 3,200 protein-class graph summaries against the preserved initial audit. All comparisons agree. The later source-linkage addendum closes the initial orientation and inner-join warnings; the old result's superseded warnings are historical, not new blockers.

| Structural universe | Exact value |
|---|---:|
| Eight-base strings | 65,536 |
| Sense strings / premature-stop strings | 59,536 / 6,000 |
| Protein-preserving classes | 3,200 |
| Classes beyond the examined reference protein | 3,199 |
| Complete / incomplete nonreference classes | 2,770 / 429 |
| Connected / intrinsically disconnected complete nonreference classes | 2,530 / 240 |
| Vertices / edges in complete nonreference classes | 48,274 / 127,587 |
| Missing sense vertices in incomplete classes | 819 |
| Observed vertices in incomplete classes | 10,411 |
| Complete classes with saturated codon main effects | 32 |
| Complete classes with nucleotide-main-effect rank below codon-main-effect rank | 701 |

Completeness means all theoretical same-protein strings within this declared eight-base library were retained in the historical strict combined-DNA >50 universe. It is not completeness of the full synonymous gene graph. Missing observations are not zero function or forbidden genotypes. The 240 disconnected complete classes are disconnected before applying any measured score, so their cross-component endpoints cannot be used as finite score-barrier examples.

`inventory.json` and `classes.csv` identify every class by its exact protein SHA-256, graph hash, membership, expected/observed vertices and edges, missing keys, component sizes, codon-factor architecture, design ranks and residual dimensions. The original hash convention uses the standard-sense translation including terminal stop; a separate biological hash corrects the fixed GTG initiator to methionine and excludes the stop. These conventions partition the library identically. Every key belongs to exactly one class, so classes share no genotype vertices. They nevertheless share one reporter, the same local mutation region, the same pooled libraries and processing. All proteins share 309 fixed residues out of 313, at least 98.7% identity at the common aligned positions. A genuine independent protein-family split is unavailable.

## Algebraic eligibility and baseline

Every theoretical class is a Cartesian product of the allowed choices in four affected codons, verified explicitly against its keys. The nucleotide graph of an individual synonymous codon factor need not be a clique. For factor sizes s_1,...,s_4 and n=product(s_j), the categorical codon-additive rank is

p = 1 + sum_j(s_j-1), with residual dimension q=n-p.

Use this codon-additive space, not the narrower nucleotide-additive space, for an inter-codon residual. The 701 rank differences show why a nucleotide residual can merely reflect an ordinary within-codon effect. Exact rational local ranks were computed; 25 representative complete design matrices independently verified the product-rank formulas. The candidate rule q>=3 concerns geometric residual dimensions, not independent biological error degrees of freedom. The 348 connected complete nonreference classes with q<3 remain inventoried but outside the proposed primary comparison; q=1 or 2 is not intrinsically meaningless, simply outside this proposed minimum-dimensional primary design.

There are **2,182 complete connected nonreference classes with q>=3**, containing 38,530 vertices. A conservative local-stratum quarantine removes the entire known-reference central amino-acid pair VP, not merely the examined protein. Six otherwise eligible classes, 192 vertices, are quarantined. The remaining **2,176 classes and 38,338 vertices** are candidates for the new analysis.

Of all 2,182 candidates, 1,721 have clique codon factors and 461 have at least one nonclique factor. After quarantine, these counts are 1,715 and 461. In a clique product, a common monotone transformation of a codon-additive field has connected nonempty superlevels. In a nonclique factor, even a codon-additive score can have a local barrier. Thus positive barriers in all these graphs cannot uniformly be interpreted as interactions between codons. Barrier diagnostics must compare the categorical codon-additive field and report the two architectures separately. No barrier implies an RNA folding mechanism.

## Proposed class grouping and sealed measurements

Group classes by their amino acids at full variable codons 157 and 158. This is a local coding stratum, not a gene family or independent experiment. Quarantine the VP stratum containing the already examined reference. Sort the remaining 399 strata by SHA-256 of `GILL_C3_V1|` followed by their two-letter amino-acid key. The first 200 strata are development, the other 199 validation. Filtering only by the structural rule above yields:

| Proposed partition | Classes | Vertices | Codon rank | Residual dimension |
|---|---:|---:|---:|---:|
| Development | 1,124 | 20,028 | 4-12 | 3-60 |
| Primary validation | 1,052 | 18,310 | 4-12 | 3-60 |
| Quarantined eligible VP classes | 6 | 192 | 8 | 24 |

This grouping keeps shared full-codon amino-acid contexts together and keeps genotype sets disjoint. It does not turn the partitions into independent families: all remain one closely related reporter collection. The primary is validation strata, not a favorable subset selected from the development result. The exact class IDs and all source hashes are in `freeze_proposal.json`, which is a proposal pending root freeze and independent review.

All additional RNA outcomes remain sealed. The previously completed reference pilot remains a failed frozen predictor, with its RNA3 outcomes already examined. It supplies neither an untouched holdout nor a positive result for this expansion. The source metadata still establish six deposited pooled libraries, not independent cultures or RNA/DNA biological pairing. CDS_1/2/3 are RNA and CDS_8/9/10 DNA. Numeric suffixes are not a pairing key. All DNA coverage information was exposed historically; this cycle reopens only sequence column A.

## Reproducibility-first endpoint

Fix c=0.5 and save the c=0.25 and c=1 sensitivity outputs without selecting between them. For each selected class form the following fields on its fixed observed genotypes:

T = [log2(RNA1+c)+log2(RNA2+c)-log2(DNA8+c)-log2(DNA9+c)]/2,

V = log2(RNA3+c)-log2(DNA10+c).

These are processed relative-library fields, not absolute concentration, protein output or fitness. Keep zeros and every declared genotype. Do not use shrunk source log-fold changes, marginal SEs or p-values. No individual DNA-positivity filter is added. The historical combined-DNA filter defines the conditional source universe.

Define the fixed uniform factorial projector P f=sum_j E[f|codon_j]-3 E[f], and Q=I-P. The projector depends only on sequence membership, not outcomes. Compute u=QT and v=QV. The primary endpoint over the 1,052 validation classes is

R = [sum_g <u_g,v_g>/n_g] / sqrt([sum_g ||u_g||^2/n_g][sum_g ||v_g||^2/n_g]).

All declared classes enter. This weights each class's mean squared contribution rather than letting larger graphs dominate solely by size; large-amplitude classes can still dominate, so report every class's signed cosine and energy plus their distributions. Do not omit negative correlations or zero-energy classes. If either total energy is zero, R is undefined and the residual test is uninformative. Individual zero-energy class cosines are undefined and retained as such, not silently removed from the primary sums.

A proposed **descriptive target R>=0.5** defines a clear pre-outcome repeatability target; it is a design convention, not a calibrated significance threshold or a universal standard. R below 0.5 fails that target; R<=0 lacks positive pooled agreement under this definition. R>=0.5 supports repeated processed-measurement organization, which can include persistent technical bias. It does not establish reproducible biological interaction. Root should adjudicate this target before any protected read; no data have been used to choose it.

The validation projection defines the residual endpoint on the held-out field. It is not a predictor fitted from validation responses. No pooled in-sample residual is evidence of predictive success. Secondary fixed, unshrunk training residual prediction u versus held-out v can be compared with the zero-residual baseline by class-weighted squared error. Report this comparison even when cosine is high; correlated patterns can have poor amplitude reproduction. Do not fit a new shrinkage parameter, structural score or higher-order model to improve it.

Secondary comparators are within-class constants and categorical codon main effects. Whole-field centering C f=f-mean(f) supplies a relative-field agreement diagnostic. Raw RNA within-class rank agreement for all three library pairs, using deterministic average ranks for ties, distinguishes a relative-abundance pattern from the residual question. DNA-only residual and rank agreement are also reported. These are diagnostics of the fixed library measurements, not estimates of biological reliability or an error-free latent target.

## Normalization diagnostic and exact dependence identity

Global additive offsets within a log-count library cancel under C and Q because both remove constants. Multiplicative normalization applied to the already pseudocounted counts is such an offset. Scaling original counts before adding a fixed pseudocount is not exactly equivalent: the pseudocount can cause genotype-dependent changes at low counts. Genotype-specific DNA effects and count errors do not cancel. A shared reference subtraction would cancel under Q, but it would not remove error in that genotype's original measurement.

With the same Q in each class, set r_i=Q log2(RNA_i+c) and d_j=Q log2(DNA_j+c). For all nine RNA/DNA fields and all 45 unordered field pairs including their nine diagonals, report the deterministic decomposition

<r_i-d_j,r_k-d_l> = <r_i,r_k> + <d_j,d_l> - <r_i,d_l> - <d_j,r_k>.

Also report its averaged-training versus validation version, using the mean of r_1,r_2 and d_8,d_9 against r_3,d_10. Save all four signed terms and their total per class, and the same aggregate with weights 1/n_g. Check numerical equality with a scale-aware floating-point tolerance. A repeated DNA denominator contributes its squared norm; different DNA libraries can also share genotype-specific bias. Compare RNA-only residual agreement with every cross-DNA scenario. The ordinal-aligned scenario is an illustrative library assignment, not verified biological matching.

This identity assumes no independence and yields no p-value. Its terms are not causal fractions, variance components or a biological noise ceiling. Nine reused-library scenarios are nine dependent fields, not nine independent replications or a confidence region. Barriers from different fields can maximize at different endpoint pairs; a range of maxima is not a fixed-pair uncertainty interval.

## Gates, failure classifications and exposure

Before any RNA read, verify source, inventory and code hashes; root must freeze the protocol and obtain independent review. Decode only the permitted selected keys/columns. Require exactly one row per fixed key and finite nonnegative integer selected counts. Invalid or missing selected rows stop execution as a source/implementation failure; they do not authorize dropping a class. Keep zero RNA counts. Excluded classes and incomplete genotypes remain outside outcome inspection. Save the complete training receipt and computations before opening any new RNA3 cells. No development outcome may change the frozen membership, primary endpoint, pseudocount, weight, transformation or benchmark list.

- Failed structural/hash/coverage gate: source, eligibility or software problem, not a biological hypothesis test.
- Zero aggregate residual energy: the intended residual agreement is unidentified/uninformative.
- R below the frozen target: failure of the declared processed-library reproducibility target, not proof against folding or a rescue of the earlier predictor.
- Dependence on DNA scenario, rank instability or inconsistent interaction summaries: a descriptive measurement discrepancy whose biological versus technical causes remain unresolved.
- R meeting the target: bounded repeatability within one reporter and these library preparations; no independence, causal mechanism, new-gene transfer or calibrated uncertainty claim.

This cycle decoded only workbook column A. It read no new RNA, DNA-count, fold-change, SE, test-statistic or p-value cells, and did not reopen the reference count vectors. The A-only reader passed a fixture containing deliberately invalid non-A shared-string indices and numeric entries, demonstrating that those cells are skipped before value decoding. The initial optional-sharedStrings handling error was corrected before the first workbook row was decoded; it did not expose outcomes. No author code, raw reads, new experiments, sequence engineering, repository/manuscript edit, commit or publication occurred.

## Reproduction and evidence

Run `work/env/Scripts/python.exe work/gill_multiclass_cycle3/inventory.py`, then `work/env/Scripts/python.exe work/gill_multiclass_cycle3/verify_inventory.py` from the workspace root. Both operate only on declared metadata and structural source material; the inventory uses column A only. The input revision is e9846dca9f74b2c5c34e984eed970e996b9df2bf. Outputs and proposal remain solely in this owned work folder.

Core evidence is `inventory.json`, `classes.csv`, `summary.json`, `manifest.json`, `verification.json`, `freeze_proposal.json` and `exposure_ledger.json`. The inventory includes every exact protein class and eligible partition, allowing an executable selected-key freeze without scanning unselected outcome cells. Existing construct/amplicon linkage, replicate metadata and the earlier failed reference protocol remain the authoritative source context. The proposed numerical target and grouping are new pre-outcome design choices, not completed biological results.
