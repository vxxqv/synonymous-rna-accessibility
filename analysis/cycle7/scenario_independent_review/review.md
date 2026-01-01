# Independent branch D preflight review

Verdict: **PASS** for the held mathematical protocol, implementation and synthetic-only design. No blocking bug or algebraic counterexample was found. No observed successor field or summary was opened. Production file-reading, source seals and execution authorization remain root's separate runner responsibility.

## Exact review point

- protocol.md: `f3e6d1046083bfc2587d96670c89a384c846429fe7747e29fffd3e663ae3cdac`
- derivation.md: `db37933a6075c830afb2ad4bf4a5aa55d7994968fcc830443a053cea8dd80e8d`
- scenario_geometry.py: `667422c85a13d5a44761b1248fa0bbae9f2469306b80f546ab9c496e3e7128cf`

The candidate-manifest file hashes are checked independently in review_manifest.json. No shared source was changed.

## Quantifier reconstruction

At a finite common absolute threshold tau, let A_r(P) mean that every vertex of path P, including its endpoints, has score at least tau in joint field r. Cycles may be removed, so finite simple-path enumeration is exact.

- Fixed route means there exists P such that every r accepts P. Its capacity is max_P min_r min_v f_r(v), exactly the capacity of the vertexwise minimum field.
- Adaptive existence means every r has some P_r. Its threshold cut is min_r max_P min_v f_r(v). It does not assert that an implementable route-selection policy knows the scenario or that routes were biologically traversed.
- Some scenario connects iff tau<=max_r C_r. Thus all scenarios are blocked iff tau>max_r C_r. A threshold above min_r C_r only establishes that at least one scenario is blocked.
- The vertexwise maximum field is an optimistic Cartesian mixing. Its capacity is at least every scenario capacity. Therefore a box obstruction implies an all-scenario obstruction; an obstruction occurring only in the box but absent from a joint scenario is impossible under this convention. The useful distinct category is the reverse: all joint scenarios blocked while the optimistic box connects.

The four capacity cuts are ordered C_fix<=C_ad<=C_max<=C_upper. E=min_r min(f_r(s),f_r(t)) is the correct endpoint eligibility cut. Restricting tau to[0,E] is essential. A scenario-specific endpoint normalization would answer a different question and can reverse interpretation.

The five half-open category intervals intersected with[0,E] are disjoint and exhaustive at each eligible threshold. Exact ties belong to the preceding category because acceptance is weak >=. Individual pairs can have nonempty intervals in several categories, so category nonemptiness counts over pairs are not a partition. Threshold widths are not prevalence or probabilities.

## Damaging controls and limits

Parallel two-route diamond with endpoints2 in both scenarios and intermediate scores(1,0) versus(0,1): C_fix=0 while C_ad=C_max=C_upper=1. A route exists in every scenario for0<tau<=1, but no one route works in all of them.

Single series path with the same two intermediate assignments: C_fix=C_ad=C_max=0 while C_upper=1. Each joint field is blocked at0<tau<=1; switching the field independently at the two vertices makes an artificial connected optimistic field. The supplied derivation states these distinctions correctly.

Disconnected base graphs have all capacities negative infinity, so every eligible nonnegative threshold, including0, is in upper_box_obstruction. Connected zero-valued graphs retain routes at0. One scenario collapses fixed/adaptive/some/upper capacities. One vertex yields zero distinct endpoint pairs rather than a fabricated self-pair. These boundaries are correctly handled.

For unknown entries, all-lower completion0 and all-upper symbolic infinity give sound coordinatewise bounds for every monotone metric. The bounded helper's conditions are sufficient universal-completion certificates; inability to certify remains unresolved. It never emits a false exact five-way record for an incomplete collection. The implementation is deliberately conservative: one unknown anywhere causes all pairs to use bounded status even if some pairs are unaffected. This may lose exact classifications but does not create false evidence. Binary/ternary finite completion tests support the implementation, while monotonicity supplies the argument for all finite nonnegative completions. Infinity is an upper bound, not an observed value.

## Independently authored checks

reference_checks.py imports the pinned implementation solely as the subject under test. The oracle independently enumerates every permutation of possible interior vertices and checks the literal path/scenario quantifiers. It does not use the production Floyd closure or production classification logic to derive expected labels.

- All64 undirected graphs on four vertices and all256 ordered pairs of binary four-vertex fields:16,384 graph/field cases.
- All98,304 pair metric records and393,216 threshold classifications agree, including exact ties0/1, interior threshold1/2 and endpoint ineligibility3/2.
- All665 incomplete two-field assignments from{0,1,None} on a three-vertex path; all unknown completions from{0,1,2}:46,683 completion-pair checks and233,415 bounded-helper checks. Defined known values never changed. Unknowns were also completed above the known range.
- A synthetic72-vertex,27-field constant-score fixture retains all2556 pairs and69,012 scenario-pair entries. The mean point is not included as a28th scenario.
- One-scenario and one-vertex controls agree with the definitions.

This finite validation does not certify every possible graph or rational landscape by enumeration. The independently reconstructed finite path argument validates the general quantifiers; the computation challenges their implementation in the stated ranges. No probability model, correlation estimate, favorable scenario pruning, biological validation or novelty claim is introduced.

## Required execution boundary

Before observed evaluation, root should freeze the candidate hashes, require all27 fields in the declared order, all72 declared sequences, exact204 edges, all2556 distinct pairs and69,012 scenario-pair records, with unknowns retained. Any source/schema/hash discrepancy, false certificate or arithmetic disagreement stops execution. This review does not authorize new outcomes itself or select a favorable subset. The mathematical preflight is ready for that root-controlled freeze.
