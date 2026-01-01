# B5 adversarial review and source-limited decision

**Verdict: correct only after an explicit scale/domain restriction on the main-marker-bound version.** Existing section 3's unrestricted-interaction impossibility and its sharp matched-marker/absolute-interaction intersection are correct on the stated Cartesian real-response domain. They supply no new empirical isolated-edit bound because a matched control/calibration contract and an absolute interaction bound are not authenticated. Retain B5 as an identification limitation, not a new fitted analysis.

This is a separately authored fresh challenge with inherited context, not a blinded external review. The existing result is `work/research_cycle8/marker_interval_review/derivation.md`, section 3. Root narrowed this task after discovering that derivation, so no competing new B5 analysis was built. Original proof, code, source records and outcomes were not modified. Only synthetic examples and pre-existing sequence metadata were inspected.

## What is the intended effect?

Fix a complete background and cell/amplicon context c. A genotype is `g=(e,m,u)`, where e is the intended allele, m the marker allele and u all remaining genomic, transcript and preparation features required by the estimand. An isolated readout effect would compare `Y_c(e,0,u)` with `Y_c(0,0,u)` while holding u fixed. Its interpretation depends on a declared response scale. A local motif and a guide label do not authenticate this intervention or background.

At ideal intensity level a local-word RNA/DNA ratio can be proportional to

`sum_(g matching word) N_gc * T_gc * V_gc * qR_gc / sum_(g matching word) N_gc * qD_gc`,

where T is RNA amount per represented genomic state, V is visibility of the chosen transcript segment and q are recovery terms. These are explanatory nuisance factors, not fitted source quantities. Library scaling and the source pseudocount scheme further distinguish a deposited ratio from this ideal intensity ratio. No complete-haplotype or isoform assignment is thereby identified.

The sequence-only mapping contains exactly 66 two-allele local states and no WT or isolated single-allele state. This was checked without reading the 132 response envelopes. The broader source audit does list one arrayed edit-only/marker-only panel, but its matched WT response and one marker label remain unresolved. It is not a matched four-state empirical validation contract for the 66 states. Physical aliases must share one outcome; reversed intended/marker labels do not create a new genotype or independent control.

## Existing proof: no-assumption failure survives the strongest simple challenge

Write a saturated chosen-scale contrast decomposition in one authenticated fixed background:

`Y_00=mu, Y_e0=mu+alpha, Y_0m=mu+beta_m, Y_em=mu+alpha+beta_m+eta_m`.

This is an identity with baseline conventions, not an assumption of additivity. Even if WT, every marker-only value and every combined value were known exactly, unrestricted eta permits `alpha -> alpha+t`, `eta_m -> eta_m-t`. The observed values remain identical. On an unbounded log scale alpha is unrestricted; on a positive raw scale a fixed positive mu at least requires `alpha>-mu`. Relative interaction restrictions, such as `eta_a=eta_b` or a bound on `eta_a-eta_b`, preserve the common shift and do not fix alpha's absolute level.

The synthetic reviewer constructs two positive, observationally identical full tables with WT and two marker-only ratios all equal to one, combined ratios 4 and 8, yet isolated edit-only ratios `1/1024` and `1024`. A third combined pair is held fixed as well. The interaction terms are assigned to **unordered physical allele pairs**, so the counterexample respects reversed intended/marker aliases rather than treating them as independent observations. Thus opposite intended-effect signs survive even perfect marker-only calibration. Unknown uncounted haplotypes add uncertainty; the impossibility already holds in the more favorable fully known background.

The fixed 66-by-92 design's rank 65 does not invalidate this argument. Rank is a property of the selected additive model, not of an unrestricted genotype-to-readout map. Under an explicitly imposed linear readout model, a proposed isolated contrast must be checked against the appropriate row space after canonicalizing physical alleles and controls. No such model is justified merely because a descriptor matrix is nearly saturated or a widened ordinal model is feasible.

## Existing sharp bound: accepted assumptions and boundary correction

For matched combined intervals `C_m=[CL_m,CU_m]`, marker-only intervals `D_m=[DL_m,DU_m]` and **absolute** bounds `|eta_m|<=I_m`, the existing real-scale result is

`alpha in intersection_m [CL_m-DU_m-I_m, CU_m-DL_m+I_m]`.

It is necessary by subtraction. For sufficiency, choose a difference d in `[CL-DU,CU-DL]` within I of alpha, then select `D=max(DL,CL-d)` and `C=D+d`. These values lie in their declared intervals and `eta=d-alpha` meets the cap. Independent choices across marker contexts establish sharpness of the intersection. This proof does not survive extra shared-control/covariance constraints without rechecking the joint possibility set. Empty intersection diagnoses incompatible assumptions/boxes, not an estimated isolated effect.

The separate known-WT version with `|beta_m|<=B_m` and `|eta_m|<=I_m` gives the corresponding interval `C_m-mu+[-B_m-I_m,B_m+I_m]` on a real/log response scale. **Required clarification:** if Y is a strictly positive raw ratio, all relevant raw cells must remain positive, including `mu+beta_m>0`, not only `mu+alpha>0`. The formula is not sufficient on that raw domain unless the marker domain is intersected with `beta_m in (-mu,infinity)` and endpoint openness is propagated.

Minimal counterexample: `mu=1`, combined interval `{1}`, `B=2`, `I=0`. The nominal interval is `[-2,2]`. Adding only edit positivity admits `alpha=2`, but then `beta=-2` and the implied marker-only response is -1. It is inadmissible. With strict positivity imposed on both single states, the sharp alpha set in this example is `(-1,1)`. The proposed repair is a domain restriction/qualification, not a change to existing empirical results. Formula (3) itself already warns about excluded zero endpoints; matched positive D intervals plus edit positivity handle this issue there.

The reviewer independently compared 2,376 single-marker alpha memberships with a finite value oracle, checked 5,184 two-marker systems and constructed 17,400 exact feasible witnesses, including both endpoints of the accepted positive example with alpha set `[1,2]`. These finite tests corroborate the algebra; they do not prove population coverage or establish a biological cap I.

## Unknown common monotone response is a different model

Under `Y_em=h(mu_z+a_e+b_m+i_em)` with an unknown common nondecreasing h, restrictions on latent `i_em` are not bounds on the observed response contrast eta. Logging the response does not generally remove h. The exact example uses latent values `(0,1,1,2)` with zero latent interaction and positive readout `h(z)=2^(z^2)` on those nonnegative values. Raw responses are `(1,2,2,16)` and base-two log responses `(0,1,1,4)`. The observed-log interaction is 2 despite zero latent interaction. This h can be extended strictly increasingly to the whole real domain, so a plateau is not needed for the counterexample.

An arbitrary common monotone h also leaves a positive latent scale gauge: scaling all latent effects by s>0 and replacing h(z) by h(z/s) leaves observed outputs unchanged. Numeric latent caps therefore need their own calibrated scale and cannot be borrowed from observed-log bounds. With a matched marker-only/combined comparison and **zero latent interaction**, a strict observed increase does force a positive latent edit effect; equality under merely nondecreasing h can conceal a nonzero effect. This limited sign implication does not identify numeric effect size, and its controls/restriction are not supplied by the current selected design.

Consequently the existing ordinal M0/M1 compatibility results and the section 3 additive observed-readout bounds must remain separate. Choosing `I=0`, a convenient B, a fitted response or an empirically favorable bound would be a new scientific assumption, not an identification result.

## Prior positioning and access limits

| Primary source | Inspected claim and relation | Limit for B5 |
|---|---|---|
| El Kassem et al. (2026), [cPRIME](https://doi.org/10.1038/s41467-026-78344-9); cached primary Methods/XML and source audit | Source already discusses marker non-neutrality and marker-confounding controls. Its endpoint is a detected local RNA/DNA motif signal. | No first marker-effect discovery. The current audit lacks matched full-state and absolute interaction calibration. Earlier correction metadata was checked on 25 October; no new correction-status claim is made here. |
| Otwinowski, McCandlish and Plotkin (2018), [Inferring the shape of global epistasis](https://www.pnas.org/content/pnas/115/32/E7550.full.pdf), PNAS 115:E7550-E7558, DOI 10.1073/pnas.1804015115 | Primary PDF pages 1-2 specify an additive latent trait followed by a monotonic nonlinear measurement function. This directly precedes the scale distinction above. | Their fitted model and error assumptions do not supply cPRIME interaction caps or calibrate an unknown h. PMC route presented a browser challenge; the ordinary publisher PDF was readable. No challenge bypass. |
| Crona (2020), [Rank orders and signed interactions in evolutionary biology](https://elifesciences.org/articles/51004), eLife 9:e51004 | Full primary article, especially “Additivity and rank orders” and “General circuits,” treats order compatibility, linear inequalities and signed interactions. | Ordinal compatibility is established prior; it does not identify an unmeasured single-state magnitude or authenticate experimental uncertainty. |
| Manski and Tamer (2002), [Inference on Regressions with Interval Data on a Regressor or Outcome](https://doi.org/10.1111/1468-0262.00294), Econometrica 70:519-546 | Publisher-indexed abstract/metadata describes identification regions under interval data and assumptions. It provides broad partial-identification context. | Direct DOI access returned 403; full theorem text was not inspected and is not invoked as a proof of the elementary B5 formula. No retry/bypass. |

Consensus was used once for bounded methodological discovery before the review. Returned causal-confounding sensitivity papers were not cited as evidence for this factorial genotype problem; assignment-sensitivity methods do not authenticate marker-interaction restrictions. The discovery receipt records that distinction. No broad survey, outcome acquisition or novel-theorem priority is claimed.

The exact increment of this review is the fresh adversarial confirmation of the existing restricted bound, a positive-raw-domain counterexample to an overbroad interpretation, and explicit separation of latent monotone and observed-log assumptions. The basic interval subtraction and shift nonidentification arguments are elementary established reasoning, not new generic partial-identification methods.

## Smallest missing record and stop decision

The most direct missing source record is a **matched WT and edit-only readout with authenticated complete genotype/background, cell/amplicon and preparation identity on the same scale**. It measures the desired chosen-scale contrast without assuming marker neutrality or zero interaction. For a true RNA/protein or isoform estimand, the corresponding product identity and visibility/recovery contract is additionally required. The existing unresolved panel cannot be silently promoted to that record.

A conditional partial-identification alternative can use authenticated matched marker-only and combined intervals plus an independently justified **absolute** interaction cap on that same scale. WT cancels from this particular difference, so a numeric WT observation is not algebraically necessary once that cap and the common-background contrast are justified. If instead only main-marker bounds are available, the known-WT/B/I version needs WT and both calibrated restrictions, with every raw-domain constraint retained. Transporting a cap from other markers or cell lines requires an explicit external transport assumption; a small observed marker-only effect does not estimate interaction.

No such absolute cap or transport calibration is authenticated in the supplied source records. Therefore there is no new empirical isolated-edit identified set to compute. Stop B5 at the proved nonidentification and conditional formula. Retain the measured local-state/context summaries as descriptive combined-state information; do not relabel them isolated intended-edit effects or causal biology.

All reviewed file hashes, synthetic receipts, source/access pointers and exposure boundaries are recorded in the companion manifest. The scope contains no new phenotype values, B1 interval reads, fits, thresholds, source sequence acquisition, manuscript edits or external actions.
