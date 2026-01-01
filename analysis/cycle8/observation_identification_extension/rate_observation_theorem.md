# Complete classification of the ten specified steady-state observations

Status: newly authored elementary derivation, self-checked with exact synthetic arithmetic; pending the requested fresh proof review. No priority claim is made. This is a fixed equilibrium observation problem, not structural identifiability from entire time trajectories.

## Domain and claim card

Fix any exactly calibrated amounts n,c,p>0 for the same RNA/protein population, compatible state and time scale. Assume the six positive finite constant parameters satisfy

    N' = s - (e+d_N)N,
    C' = eN - d_C C,
    P' = kC - d_P P,

at equilibrium (n,c,p). All relevant dilution must already be included consistently in the effective losses. No feedback, reverse import, maturation delay, hidden inactive compartment, heterogeneous kinetic population or variable observation gain is included. k is completed protein production per cytosolic RNA per time, not an unqualified initiation rate.

The candidate additional scalar observations are exactly

    s, a=e+d_N, e, d_N, d_C, J=en, k, d_P, H=kc, f=e/a.

J is RNA export flux. H is completed protein production flux. They are different quantities. An observation is an exact independently authenticated scalar with its stated calibration, not an annotation saying a similarly named assay was run. Repeated/derived scalars can be redundant. This classification makes no claim about every possible nonlinear or discontinuous measurement outside this list.

For any subset S of these ten labels, call it globally sufficient if its observation map is injective on the entire positive equilibrium fiber. Equivalently, every admissible exactly consistent observed vector has a unique positive six-parameter solution. Empty/inconsistent observation fibers are rejected, not called identified.

**Theorem.** S is globally sufficient if and only if it contains at least one observation from B and observations from at least two distinct nuclear classes among A,E,D,F below. Its inclusion-minimal sufficient sets consist of exactly three scalars: one representative from B and one representative from each of two distinct nuclear classes. There are exactly 51 such sets, 805 sufficient subsets and 219 insufficient subsets among all 1,024 subsets.

## Normalized equivalence classes

Write x=e, y=d_N and z=d_P. The complete equilibrium fiber is

    x,y,z>0;
    s=n(x+y), d_C=nx/c, k=pz/c.

This follows in both directions from the three zero-derivative equations. The change of variables is bijective, so identifying these three coordinates identifies all six parameters.

| Class | Candidate labels | Normalize each measured value to | Class size |
| --- | --- | --- | ---: |
| A | s,a | s/n=a=x+y | 2 |
| E | e,d_C,J | e=(c/n)d_C=J/n=x | 3 |
| D | d_N | d_N=y | 1 |
| F | f | f=x/(x+y), strictly between 0 and 1 | 1 |
| B | k,d_P,H | (c/p)k=d_P=H/p=z | 3 |

These equivalences use the already known calibrated n,c,p and the equilibrium model. They are not empirical assertions that arbitrary measurements of those quantities agree. If multiple representatives of one class are supplied, all normalized values must be exactly equal. With noisy observations these become consistency/uncertainty questions outside the exact theorem.

## Global reconstruction for every nuclear pair

Let A,E,D,F denote the normalized observed values of the corresponding classes. Values A,E,D,B must be positive, while 0<F<1. The following six inverses cover every pair of distinct nuclear classes.

| Nuclear pair | Reconstruct x=e | Reconstruct y=d_N | Additional consistency condition |
| --- | --- | --- | --- |
| A,E | E | A-E | A>E |
| A,D | A-D | D | A>D |
| A,F | AF | A(1-F) | None beyond stated domain |
| E,D | E | D | None beyond stated domain |
| E,F | E | E(1-F)/F | None beyond stated domain |
| D,F | DF/(1-F) | D | None beyond stated domain |

Any B representative then gives z. Recover s=n(x+y), d_C=nx/c and k=pz/c. Every division is by a positive known quantity or by F or1-F, both nonzero in the domain. The reconstructed rates are positive precisely under the table's conditions. Direct substitution proves the six inverses globally, without a local rank argument. For a larger S, reconstruct from any pair and B representative, then check every other supplied observation. Consistent values yield the same solution; inconsistent values yield no solution.

## Necessity and minimality

If S contains no B representative, changing z to any other positive number leaves all selected observations fixed but changes k and d_P. Therefore it cannot identify all rates, even if every nuclear observation is known.

If S includes B but has no nuclear class, x and y are free. If it has only one nuclear class, a positive one-parameter ambiguity remains at every point:

| Only observed nuclear class | Positive indistinguishable alternatives |
| --- | --- |
| A | Change x inside(0,A), set y=A-x |
| E | Keep x, vary y>0 |
| D | Keep y, vary x>0 |
| F | Scale both x and y by any common positive factor |

Additional representatives of that same class impose only its redundancy equalities. Thus the criterion is necessary at every admissible positive point, not merely generically. Together with the inverse table it proves the theorem.

Every sufficient set contains a B representative and a pair of distinct nuclear representatives; removing any one of those three destroys sufficiency. Conversely a larger sufficient set contains such a triple and is not inclusion-minimal. Nuclear representative pairs number

    2*3 + 2*1 + 2*1 + 3*1 + 3*1 + 1*1 =17.

Multiplication by the three B choices gives 51 minimal triples. There are2^7=128 nuclear subsets. Exactly1+(2^2-1)+(2^3-1)+(2^1-1)+(2^1-1)=13 occupy zero or one nuclear class, leaving 115 sufficient nuclear subsets. The nonempty B subsets number 7, giving 115*7=805 sufficient subsets and 1,024-805=219 insufficient ones. `observation_subsets.json` lists every subset and all 51 triples, making the finite combinatorial enumeration inspectable. This count is a consequence of the proof, not its replacement.

## Boundary, calibration and practical stability

The theorem excludes d_N=0 and f=1. Those can be meaningful alternative model assumptions, but must not be inferred from a nearly cancelling fitted quantity or borrowed from another study. Under an assumption d_N=0 known in advance, one A/E representative and one B representative suffice; that is a different information problem. An observed A=E is inconsistent with the present strict-positive model. Zero amounts, zero protein loss/production and unequal measurement scales can create other boundary models and are not silently included.

Global uniqueness under exact inputs does not guarantee practical stability. For A,E, y=A-E is a near-cancellation when y is small: perturbations satisfy |Delta y|<=|Delta A|+|Delta E|, so relative error can be large compared with y. For A,F and E,F, the inferred y includes1-F; absolute fraction error is amplified relative to y as F approaches1. For D,F, x=DF/(1-F), so relative fraction sensitivity is1/[F(1-F)]. These identify an actual conditioning issue, not a reason to discard an unfavorable measurement.

Not every sufficient set is intrinsically ill-conditioned at d_N approaching0. E,D directly supplies x,y; A,D obtains the non-small x=A-D near that boundary. Their practical performance still depends on the measurement precision of a small D, but they avoid the same subtraction or division by 1-F. The opposite boundary e approaching0 has its corresponding sensitivities. Uncertain n,c,p and cross-assay calibration add further errors not represented by this exact exercise.

## Occupancy with unknown residence time

Suppose, additionally and only hypothetically, a calibrated occupancy O obeys the ideal steady-flow relation O=H*tau for one unknown residence time tau>0, with no unaccounted abortive pool. Per-RNA occupancy rho=k*tau is equivalent when c is known, since O=c*rho. These equalities are not established by a heavy-polysome metafraction score.

For every z>0, H=pz and tau=O/(pz)>0 give the same positive O. Consequently this added occupancy does not replace any B observation or reduce the rate-identifying criterion. If S already identifies all six rates, O also identifies tau. For identification of the augmented rates-plus-tau parameter vector, O plus a sufficient original S is necessary and sufficient; its minimal sets are the same51 triples with O appended. Without O, tau is unobserved. Known tau together with O would supply H, but a borrowed residence time is an additional assumption. O=0 is inconsistent with the strict-positive ideal relation, not evidence for a zero completion rate in this model.

The statement concerns the stipulated relation O=H*tau. It does not equate the RNA export flux J with protein production H, infer initiation from completion, or apply a conservation identity to stalled/aborted occupancy without checking its assumptions.

## Verification and novelty limit

The exact checker enumerates all 1,024 subsets. At five positive rational anchors, including near each nuclear boundary, every sufficient subset reconstructs all rates and every insufficient subset receives an explicit positive alternative with identical selected observations. Boundary/redundancy failures are tested. The code's first passed draft was preserved before adding explicit domain refusal for F outside(0,1); no theorem, model or source-dependent selection changed. No rate or biological observation was fitted.

Müller et al.2024 already gives the nuclear/cytosolic subsystem and obtains extra information from compartment-resolved labeling trajectories. The protein balance and algebraic elimination are standard. Chis et al.2011 distinguishes global uniqueness from local identifiability. The table is a complete elementary classification for this specified candidate list and observation model, with no novelty-priority claim beyond inspected sources. Fresh adversarial proof review remains required before manuscript promotion.
