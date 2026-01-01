# Observation model and recovery-identifiability claims

Status: internally derived here, pending separate proof review. No observed outcome was read or fitted. Claims are elementary finite-factorization/interval arguments, with no mathematical-priority assertion.

## Assumptions and objects

Let V be a finite set of genotypes with reference w. Let r index output preparations and t index technical input preparations. D_v>0 is source plasmid abundance. Let e_vr>0 describe effective delivered-template representation relative to D_v, a_vr>0 total relevant RNA per effective template, p_vr in[0,1] the fraction of that RNA carrying the fully included target exon, and q^O_vr>0 output RT/PCR/sequencing recovery. Let q^I_vt>0 be technical input recovery. A literal inclusion/(inclusion+skipping) PSI equals p only if the declared RNA population consists of those two outcomes or all other isoforms are properly accounted for.

At the intensity/probability level, assume

`mu^I_vt = alpha_t D_v q^I_vt`,

`mu^O_vr = beta_r D_v e_vr a_vr p_vr q^O_vr`.

The nuisance scales alpha,beta may include sequencing depth and shared preparation yield. Here p is a fraction in the explicitly defined RNA population, not a kinetic splice-choice rate; differential isoform lifetimes can separate those quantities. No Poisson/multinomial model or calibrated sampling uncertainty is inferred merely by writing intensities. In particular, an observed ratio of random counts is not an unbiased estimate of the ratio of their means in general.

Require p_wr>0 and positive reference/input intensities whenever a ratio is used. Then ideal reference-relative input-normalized included-RNA enrichment is

`R_vrt = (mu^O_vr / mu^O_wr) (mu^I_wt / mu^I_vt)`

`= (p_vr/p_wr) B_vrt`,

where

`B_vrt = (e_vr a_vr q^O_vr)/(e_wr a_wr q^O_wr) * q^I_wt/q^I_vt`.

Here B_wrt=1. Source plasmid abundance and common per-preparation scales cancel. Effective intracellular representation, total RNA abundance/stability and differential recovery do not generally cancel. The count-level full-library totals also cancel algebraically in the separately normalized relative ratio. Historical median-input normalization is a different estimator; its source reconstruction and missingness rules remain governed by the frozen FAS protocol/amendment.

## Claim 1: nonidentification under unrestricted recovery

For one positive finite reference-relative field y with y_w=1, any constant inclusion p_v=p0 in(0,1] can reproduce y: set e=a=q^I=1, q^O_w=1 and q^O_v=y_v. If primitive recoveries are required to be probabilities at most1, instead take q^O_v=c*y_v with0<c<=1/max(y); their relative ratios are unchanged. Thus the same ideal observed ordering can be produced by an arbitrary positive persistent genotype recovery factor with unchanged true inclusion. No independent normalization of the same input/output product identifies which factor changed.

More generally, for any desired positive relative-inclusion field f with f_w=1, define B_v=y_v/f_v. Choose p_w sufficiently small that p_v=p_w f_v<=1 for every finite f_v. Then the same y is compatible with f. Hence an arbitrary positive ordering neither identifies inclusion effects nor their interactions without constraints on B. This is a parameter-identifiability result, not a statement that the actual assay's bias is large or that its data are erroneous.

**Preparation qualification.** One persistent B_v and a preparation-invariant inclusion field imply the same ideal R_v in every preparation after common scales cancel. They do not reproduce an arbitrary matrix of different ideal preparation fields. For example, normalized scores(1,2) and(1,1/2) for reference/variant cannot both equal one persistent B_v times one fixed constant-inclusion field. Preparation-specific B_vr could reproduce both, but that is a different model. Sampling fluctuations are also a different layer and require an explicit stochastic model. Do not treat27 observed cross-products as27 freely and independently adjustable ideal fields.

## Claim 2: exact interpretation of an anchored multiplicative expansion

Let `B=product_v[L_v,U_v]` be a **formal** nonnegative empirical box, with finite L, possibly infinite U, reference interval[1,1], and only finite candidate field values. Let M be the named nonnegative codon-unary/common-nondecreasing-response class, and `M1={f in M:f_w=1}`. For kappa>=1, define

`Q_WT(kappa)={b>0:b_w=1, kappa^-1<=b_v<=kappa}`.

Define an anchored expansion with interval[1,1] at w and `[L_v/kappa,kappa U_v]` elsewhere. Then

`exists f in M1, b in Q_WT(kappa), b*f in B`

is equivalent to M1 intersecting that anchored expansion.

**Proof.** Necessity follows by dividing each coordinate by b. For sufficiency, for each f_v>0 choose a value y_v in the nonempty intersection `[L_v,U_v] intersect [f_v/kappa,kappa f_v]`, then set b_v=y_v/f_v. At f_v=0, expanded feasibility forces L_v=0; choose y_v=0 and any allowed positive b_v. At w choose y_w=f_w=b_w=1. Bounds at infinite U are interpreted without division by infinity. The chosen b is one persistent factor per genotype for this one chosen field, not one per preparation. QED.

The conditional minimum kappa is a mathematical model-rescuing WT-relative allowance on the combined inclusion-extrinsic multiplier. It is not an estimate, confidence interval or authenticated bound on actual recovery. This multiplier includes effective template and total RNA effects, which may be genuine biology rather than technical errors. Bounding instrument/RT/PCR recovery alone is insufficient if genotype-dependent total RNA remains unrestricted. Treating B as containing the true intensity field is an additional assumption unsupported by an empirical range alone.

## Claim 3: all-coordinate expansion instead measures pairwise spread, with a square

Assume in addition that M is closed under positive scaling: `f in M => c f in M` for c>0. The common-response codon-unary class has this property because c*h is nondecreasing. Let

`B_lambda=product_v[L_v/lambda,lambda U_v]`, including the reference,

and

`Q_pair(K)={b>0:b_w=1, max_v b_v/min_v b_v<=K}`.

For every finite lambda>=1,

`M intersects B_lambda`

if and only if

`exists f in M1, b in Q_pair(lambda^2), b*f in B`.

**Proof, forward.** Choose f' in M intersect B_lambda. For every positive f'_v, select y_v in `[L_v,U_v] intersect [f'_v/lambda,lambda f'_v]`, and put d_v=y_v/f'_v in[1/lambda,lambda]. At f'_v=0 choose y_v=0,d_v=1. The reference has f'_w>0 and y_w=1. Set f=f'/f'_w and b=f'_w d. Then f_w=b_w=1, b*f=y, and the ratio max b/min b equals max d/min d, at most lambda^2.

**Proof, reverse.** Given f,b with pairwise spread at most lambda^2, let `c=sqrt((max b)(min b))`. Then d=b/c lies in[1/lambda,lambda]. Let f'=c f, which belongs to M. Since y=b*f=d*f' lies in B, f' lies in B_lambda. QED.

Consequently the **infimal** pairwise spread allowance equals the square of the infimal all-coordinate symmetric expansion factor, under these assumptions. This does not assert attainment, rationality or statistical calibration of either infimum. An explicit feasible f' also yields a persistent-factor construction by coordinatewise clamping into the original box and the normalization above. That construction fits a field in the rectangular box; it need not equal any observed joint scenario.

Lambda is a factor at least1, not a percentage perturbation. If a reported radius is delta with lambda=1+delta, the corresponding pairwise factor is(1+delta)^2 and its excess above1 is2*delta+delta^2. For a log-radius eta=log(lambda), the pairwise log-spread is at most2*eta. None is a measured biological recovery percentage.

Two bound conventions must not be interchanged: a WT-relative bound b_v in[1/kappa,kappa] permits pairwise spread up to kappa^2. A pairwise spread cap K is stronger than that WT-relative cap with kappa=K. If primitive unnormalized recovery factors are each within a factor c of a common baseline, ratios to WT may lie in[c^-2,c^2]. These squares are normalization consequences, not additional biological evidence.

For the two corresponding expansion optima, scaling closure gives `lambda_all <= kappa_WT <= lambda_all^2`, interpreted as bounds on infima. The first inclusion follows because the anchored expansion is contained in the all-coordinate expansion; the second follows from the pairwise-spread construction and b_w=1. A toy constant-field model with original box point(1,K) has lambda_all=sqrt(K), while the anchored WT-relative allowance is K, showing why the two metrics can differ. This toy example is not a claim about the actual FAS optimum.

## Claim 4: a rectangular completion is not simultaneous persistent correction

Claim 2/3 contain quantifiers `exists ONE y in B, ONE f, ONE b`. They do not establish either

`exists b,f such that y^s=b*f for ALL observed scenarios s`,

or the weaker but different condition

`exists b such that y^s/b belongs to M for ALL s`.

They also do not follow by fitting each scenario with its own b^s. The two positive fields in Claim 1's preparation example are separately corrected to constant inclusion by different factors2 and1/2, but no one factor makes both constant. This is a minimal counterexample for the constant-inclusion explanation, not a proof that the full FAS codon-unary class cannot fit those fields.

An empirical box discards cross-preparation dependence and includes mixtures of extrema never observed together. To interpret an expanded-box optimum as persistent bias explaining every preparation, an additional model must link the ideal field, preparation noise/interactions and observed scenario vectors. No such model is authenticated by the existing27 empirical cross-products.

## Claim 5: additive radii are not fold-recovery bounds

An additive expansion `max(0,L_v-epsilon)<=f_v<=U_v+epsilon` concerns absolute score-unit differences. It does not impose a uniform bound on y_v/f_v. If f_v>=m>0 and |y_v-f_v|<=epsilon<m, then y_v/f_v lies in[1-epsilon/m,1+epsilon/m], so a WT-relative symmetric factor `max(1+epsilon/m,1/(1-epsilon/m))` suffices. Without a positive lower bound, the multiplicative distortion can be arbitrarily large. Conversely a fold factor kappa implies absolute differences at most(kappa-1)M only if f_v<=M is known.

For an anchored physical relative-score error, WT must remain1. An all-coordinate additive expansion permits the corrected reference to move; renormalizing it changes absolute errors at other coordinates. Therefore its formal radius cannot automatically be called an anchored additive recovery allowance. The exact anchor/gauge convention must accompany the number.

A positive finite recovery multiplier preserves exact zeros. Therefore multiplicative widening leaves a zero interval[0,0] unchanged, whereas additive widening can allow positive f. Neither operation turns a sampled count zero into evidence of true zero inclusion. Unknown[0,infinity] coordinates remain noninformative under finite widening. Any biological interpretation of an additive or multiplicative radius needs these distinctions and an observation-noise assumption stated explicitly.

## Dependency and evidence labels

Source-supported leaves: pooled included-band assay, technical inputs and repeated transfections, and incomplete direct-PSI overlap, from authenticated 2016/2019 Methods and earlier gates. Unchecked biological leaves: proportional effective template sampling, constant totalRNA/recovery across variants, persistent versus preparation-dependent nuisance, two-isoform completeness and true-field coverage by the empirical box. Internal mathematical leaves: factor cancellation, coordinatewise interval intersection and positive-scale closure. Synthetic examples corroborate boundary/counterexample arithmetic but do not prove the universal claims. No new molecular mechanism or foundational theorem is asserted.
