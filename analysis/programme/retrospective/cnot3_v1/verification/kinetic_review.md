# Independent challenge of the two-rate identified set

**Primary verdict: proved as written for the declared exact model and finite closed rectangular admissible sets.** The author derivation correctly distinguishes the basic algebraic model from the model requiring a nonnegative perturbed background. Its sharpness argument is valid. No correction to the held author sources is required.

The accompanying claim that two paired rates do not generally point-identify a direct pathway contribution is supported by explicit distinct-parameter counterexamples. This is not a claim that every boundary instance is nonidentified: k0=0 forces p=0, and additional externally stipulated restrictions can collapse an interval to a point. Biological use remains conditional on the kinetic endpoint and decomposition assumptions.

## Reviewed claim and provenance

For fixed nonnegative k0,k1, bounds 0<=rL<=rU<=1 and dL<=dU, the model is k0=b+p and k1=b+delta+r*p, with b,p>=0. Unknown r and delta range over their whole declared intervals. The physical variant additionally imposes b+delta>=0. The conclusion is the exact projection of this feasible parameter set onto p, with constructive attainment of every retained value.

I first reconstructed the proof and a distinct exact calculation without opening the author lane. The completed derivation, code and results were hashed in `precomparison_freeze.json`. Only then did I read `kinetic_model/derivation.md`, `kinetic_bounds.py` and its manifest. This is a separately authored review with inherited project context, not a blinded experiment or a claim of model-level statistical independence.

Audited author hashes:

* `derivation.md`: `c9e2a82f62c97949b01f87d541ccdf0ed6ca1e1445f5e9f0cb288471b85cc2b9`
* `kinetic_bounds.py`: `1b5f4b33f40d71dc5331bdf367bbb9839dbd14000ece0d633df0e04983899664`
* `verify.py`: `7d89361c6f9aac5944ba40e0ce08f0121c955d8b0583964ef657219b65a0bc80`

Their hashes and the additional author receipt hashes were unchanged at completion of the comparison.

## Independent critical argument

Put D=k0-k1, aL=1-rU and aU=1-rL. For fixed p, eliminating b leaves q=(1-r)p=D+delta. Basic feasibility is exactly

    0<=p<=k0,
    [aL*p,aU*p] intersects [D+dL,D+dU].

The physical condition is q>=p-k1. Therefore the complete sets are

    P_basic = {p: 0<=p<=k0, aL*p<=D+dU, aU*p>=D+dL},
    P_physical = {p in P_basic: rL*p<=k1, p<=k0+dU}.

Choosing q from the indicated intersection gives r=1-q/p, delta=q-D and b=k0-p when p>0. At p=0, the intersection forces q=0 and delta=-D, with arbitrary allowed r. This reconstructs a witness for every point and proves sharpness without dividing by a coefficient that might vanish. It is identical in content to the independently completed author argument.

For rectangular rate uncertainty I additionally checked the author's elimination by a different representation. Let c=r*p and keep (p,b,c). The first rate requires L0<=b+p<=U0. The second rate can be met by some delta in [dL,dU] exactly when

    L1-dU <= b+c <= U1-dL.

With 0<=p,b and rL*p<=c<=rU*p, these are a bounded polytope. Physical nonnegativity adds b+dU>=0 and c<=U1. Projection onto p agrees with the author's formula. The equivalent interval-sum proof in the author report is sound, including the necessary check that the physically truncated background interval is nonempty.

## Obligation matrix

| Obligation | Verdict | Evidence or qualification |
|---|---|---|
| Signs and model identity | Passed | D=(1-r)p-delta follows directly. Positive delta reduces observed net suppression; negative delta can produce a positive D even at p=0. |
| Exact-rate interval | Passed | Independent interval-intersection proof above. |
| Physical constraints | Passed | Both rL*p<=k1 and p<=k0+dU are required in the bounded model and are included. |
| Sharpness | Passed | Explicit witnesses cover every retained point, not only interval endpoints. |
| r=1 and zero coefficients | Passed | Feasibility requires delta=-D; basic [0,k0], physical [0,min(k0,k1)]. Code treats zero coefficients as feasibility conditions. |
| Negative D, empty sets and zero baseline | Passed | No clamping; contradictory restrictions remain empty. Compensation can allow D<0 with positive p. k0=0 forces p=0. |
| Finite rectangular rate bounds | Passed | Independent 3D polytope representation plus interval argument and exact checks. |
| Paired-rate and correlated uncertainty | Passed with stated scope | The author correctly limits sharpness to the declared rectangle. Correlated subsets inherit only outer bounds unless their actual joint constraints are retained. |
| Half-life limits and units | Passed with stated scope | Reciprocal bounds reverse endpoints. If all lifetimes must remain finite, a zero rate from an infinite upper bound is an infimum, not an attained value. An unknown common time scale prevents absolute calibration. |
| Compensation interpretation | Passed | The author distinguishes suppressed contribution S, baseline contribution p and net response D, and does not identify them without restrictions. |
| Causal pathway attribution | Conditional | Requires valid exponential/comparable endpoints, residual-activity restrictions and bounds on other changes. No biological data validate these assumptions here. |
| Cross-gene common parameters | Out of scope of the single-pair sharpness statement | If one common r or another shared latent quantity is required across transcripts, per-transcript intervals need not be jointly attainable. Retain this qualification in any later application. |

Dependency chain: declared model and allowed parameter set -> elimination of b -> intersection of q intervals -> physical half-line -> constructive parameter recovery -> endpoint interpretation. Every mathematical leaf here is elementary algebra. No external theorem or finite enumeration is needed for the universal proof.

## Decisive counterexamples and what they establish

1. **Nonidentification:** k0=1,k1=1/2,delta=0 admit (b,p,r)=(1/2,1/2,0) and (0,1,1/2). Both backgrounds are physical. The same observations have different baseline pathway contributions.
2. **A physically impossible basic fit:** k0=1,k1=0,r=1/2,delta=-1/2 force p=1,b=0. The basic fit exists but the perturbed background is -1/2, so the physical model is empty. This challenges the omission of the physical constraint, not the reviewed author result, which includes it.
3. **Dependence cannot be replaced by marginal rate differences:** if the true allowed pairs are (x,x), 1<=x<=2, with r=0 and delta=0, p is exactly zero. Independent marginal rate intervals admit (2,1), hence p=1. The rectangular outer approximation can therefore be substantially wider than the true identified set.

Useful restrictions are correspondingly precise. Known r<1 and known delta can identify p, provided every constraint holds. With bounded compensation and residual activity, the sharp interval may provide information, but it need not. Bounds are only scientifically informative when their input restrictions are externally defensible; choosing them after inspecting outcomes would not validate the mechanism.

## Executable finite evidence

`exact_reference.py` compared the independent scalar interval with a distinct 2D exact-rational polygon representation using z=r*p. All 2,880 cases passed: 1,598 nonempty and 1,282 empty. Original-equation witnesses were checked at each nonempty interval's endpoints and midpoint. The tested rate grid is {0,1/2,1,2}; residual interval endpoints are in {0,1/2,1}; shift interval endpoints are in {-2,-1,0,1,2}; all ordered bounds and both model readings are retained.

After the independent freeze, `compare_author.py` checked the held author implementation against those 2,880 exact-rate answers and 162 additional rectangular-rate cases from an independent 3D vertex oracle. All passed. It also checked 478 author-constructed witnesses directly against the original equations and bounds. No biological outcome or author result file was imported as an oracle. The author's manifest and recorded result bytes were read only to authenticate the review point.

Run from the workspace with the existing Python runtime:

    python -B work/decay_cycle4/kinetic_independent_review/exact_reference.py
    python -B work/decay_cycle4/kinetic_independent_review/compare_author.py

These finite results corroborate the proof over the stated grid only. They are not proof of all parameter values and do not establish biological validity. The output receipts are `independent_results.json` and `comparison_results.json`.

## Remaining limitations and stopping decision

There is no unresolved mathematical gap in this bounded model. The strongest safe claim is a sharp identified set conditional on explicit activity, compensation, measurement and physical-background assumptions. It is ordinary competing-rate and partial-identification algebra, not a new generic theorem. A useful next scientific step lies in authenticating measurement semantics and defensible input bounds, not in extending this algebra. No further theory, fitting or data acquisition is warranted in this review lane.
