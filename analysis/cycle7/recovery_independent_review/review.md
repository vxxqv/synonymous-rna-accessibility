# Independent recovery-identifiability proof review

Primary verdict: **proved as written under the stated model assumptions**, with one minor toy-example parameter clarification. No material algebraic or quantifier defect was found in the reviewed proof. This is an internal mathematical review, not biological calibration or validation of the observation model.

Reviewed source: claim_and_proof.md SHA256 `18c6c963a14eb70aae94d09e12335ec05bced491abae86faf313e81b51b7e771`. No observed successor field or result was opened. Source Methods, new outcomes and author synthetic functions were not re-read/imported for these checks.

## Defect/clarification first

Claim3's toy constant-field point(1,K), with lambda_all=sqrt(K) and kappa_WT=K, should explicitly say K>=1 (preferably K>1 to illustrate a difference). For0<K<1 those expressions require replacing K by1/K. The preceding use of K as a spread cap naturally implies K>=1, so this is a minor precision issue rather than a defect in the main equivalence.

## Claim1: constant-inclusion nonidentification

The intensity factorization is algebraically correct: source plasmid D_v and common preparation scales cancel, but effective template representation e, total RNA a and relative output/input recovery remain. With constant p>0 and e=a=qI=1, any one finite positive field y normalized by y_w=1 is produced by relative qO=y. If primitive qO is probability-bounded, multiply every qO by one positive c<=1/max y. This leaves all relative ratios unchanged and ensures0<qO<=1. Thus probability bounds alone do not repair identifiability. For an arbitrary desired finite positive relative-inclusion field f, choosing p_w<=1/max f and B=y/f is also valid.

This is an existence of alternative parameters, not evidence of a large bias or error in the assay. The preparation qualification is necessary and is present: one fixed f and one persistent b produce one ideal field, not different ideal fields across preparations. Differences in sampled counts are not themselves proof that the ideal field changed, because a sampling layer has not been specified. For a full r-by-t matrix, output and input factors have structure; the text correctly forbids treating the27 cross-products as freely independently adjustable ideal fields.

## Claim2: anchored expansion exactly matches one-completion bounded bias

Fix finite kappa>=1. For each coordinate with f_v>0, existence of a positive b_v in[1/kappa,kappa] and y_v=b_v f_v in[L_v,U_v] is equivalent to

`[L_v,U_v] intersect [f_v/kappa,kappa*f_v]` being nonempty,

which is equivalent to `L_v/kappa<=f_v<=kappa*U_v`. This follows directly from the two cross-interval endpoint inequalities. An infinite upper bound creates no difficulty because a finite y can be chosen. If f_v=0, either condition holds exactly when L_v=0; y_v=0 and any permitted positive b_v work. At w the fixed interval[1,1] and normalized class force y_w=f_w=b_w=1. Taking products across the finite coordinates proves both directions, with one b per genotype for the selected field.

The result needs valid nonnegative intervals, a fixed class M and a finite candidate field. It does not require any graph, probabilistic independence, paired covariance or distributional model. It does require the extra physical assumption that some true intensity field belongs to the empirical box before it can be used as a bound on an actual multiplier. That assumption is explicitly not authenticated here.

## Exact counterexample to a simultaneous-correction interpretation

Use the full binary codon Q2 with feature rows X00=(0,0), X01=(0,1), X10=(1,0), X11=(1,1), reference00. Consider positive normalized fields

`y_plus=(1,2,2,1)` and `y_minus=(1,1/2,1/2,1)`.

Their Cartesian envelope contains the constant field1, so it admits the common-monotone codon-additive model with anchored allowance1. But allowance1 forces b=1, and y_plus itself cannot belong to that class: its strict orders00<01 and11<10 require positive latent differences(0,1) and(0,-1), whose sum is zero. Hence no bounded persistent b at this allowance makes both fields model-compatible, even if separate common monotone response functions are allowed per field. This is a counterexample for the full class, not only for a constant-inclusion submodel.

More generally, for any specified finite WT-relative cap kappa, replace2 by alpha>kappa^2 and1/2 by1/alpha. After any correction within that cap, both high intermediate values exceed kappa, while the corrected endpoint11 is at most kappa and endpoint00 is1. The reciprocal-order contradiction remains. The box still contains the constant1 field. This does not claim impossibility with an unrestricted unbounded b; the cap is essential to the counterexample.

Therefore Claim2's existential one-completion statement cannot be promoted to either one b and one f explaining all fields, or one bounded b making every field separately model-compatible. Claims1/4 already distinguish those quantifiers correctly. Merely fitting fields with their own b_r also does not establish one persistent correction.

## Claim3: all-coordinate allowance and pairwise spread

Forward construction is valid, including zero coordinates. Since the reference interval is[1,1], a finite lambda gives positive f'_w. Choose d_v in[1/lambda,lambda] with y_v=d_v f'_v in the original interval, then normalize f=f'/f'_w and b=f'_w*d. At w, y_w=1 implies d_w=1/f'_w, so b_w=1 exactly. Scaling closure retains f in M. The factor spread is at mostlambda^2.

Conversely, if normalized positive b has max b/min b<=lambda^2, choosing c=sqrt(max b*min b) puts d=b/c in[1/lambda,lambda]. Scale f to f'=c f and recover membership in the all-coordinate expansion. Thus the finite-allowance equivalence and squared relation between infima hold. Positive scaling closure is essential; it is satisfied by arbitrary common nondecreasing responses. The proof does not require either infimum attained or rational gamma, and it does not turn a spread into a calibrated error percentage.

The inequalities lambda_all<=kappa_WT<=lambda_all^2 follow from set inclusion and b_w=1. WT-relative cap kappa and pairwise-spread cap kappa are different: the former permits spread up tokappa^2, while the latter implies each WT-relative factor lies in[1/kappa,kappa]. Primitive-factor normalization introduces another potential square exactly as stated.

## Claim5 and safe radius label

The additive-to-multiplicative comparison requires f bounded away from zero; without it, a fixed absolute error gives arbitrarily large relative distortions. A multiplier preserves exact zeros, so a sampled count zero cannot be treated as a proven biological zero. Renormalizing a moved reference changes additive errors. These qualifications are correct.

Recommended label for the impending quantity:

**Minimum WT-anchored multiplicative expansion that permits one common-monotone codon-additive field in the Cartesian empirical interval box.**

Equivalently, under the formal observation model, it is the smallest permitted WT-relative combined inclusion-extrinsic multiplier capable of explaining one completion of that box. It is not estimated actual recovery, a confidence interval, simultaneous persistent correction of all27 preparations, pure RT/PCR bias, a causal interaction strength, biological function improvement or validation of a source-calibration model. Effects on total RNA/effective delivery may be real biology and belong to the combined multiplier unless independently measured.

## Synthetic evidence and scope

checks.py uses separately authored exact arithmetic only. It verifies80 one-coordinate interval/bias equivalences including zero boundaries, probability-bounded constant-inclusion factorization, the Q2 cancellation counterexample, a lambda3 versus anchored9 square example, and additive half-unit errors with multiplicative factors6,51,501 near zero. These finite examples corroborate the algebra but do not prove its universal form. The proof reconstruction above supplies the general argument.

No proof was propagated into the manuscript, no scientific model fitted, no observed threshold computed, and no shared source edited. The remaining task is faithful label/quantifier transfer during root's frozen numerical and manuscript integration.
