# Cycle 7 finite-scenario accessibility protocol, candidate 1

Status: pre-outcome design and synthetic verification only. Root approval and an immutable freeze are required before any new observed derived comparison. Existing Cycle 6 exclusions, endpoints and negative conservative-barrier result remain unchanged.

## Question and fixed scope

On the completed 72-vertex, 204-edge native included-protein class, distinguish existence of one route accepted under every joint preparation field from existence of a potentially different route in every field. Also distinguish obstruction in every deposited joint field from obstruction in the larger box that permits vertexwise mixing of preparations. This is a retrospective descriptive analysis of exposed measurements. It is not independent biological validation, a causal attribution, or a new test of the global FAS exclusion.

Use all 27 BR-by-TR fields in the fixed order BR1/TR1 through BR9/TR3 and all 2,556 unordered distinct endpoint pairs. The mean corrected source-formula point is not a 28th preparation scenario and is excluded from this branch by definition. No field, endpoint, threshold or subset may be selected after inspecting the derived comparisons. Source file: `results/cycle6/full_result_v1.json.gz`; exact decompressed SHA-256 `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466`; compressed SHA-256 `e3c3d061052643046de58dc789866899b0fb14ade83f06170e2cf32e9e54e145`. Sequence contract SHA-256 `f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38`.

Graph edges are the fixed Hamming-one edges on those 72 sequences. A path includes its endpoints. Acceptance is score greater than or equal to a common absolute threshold. Scores and calculations are exact nonnegative rationals. Zero is a valid score. There is no tolerance, jitter, fitted correlation, probability or independence model.

## Primary quantities

For endpoints s,t and field r, let C_r be its widest-path capacity. Set E=min_r min(f_r(s),f_r(t)); C_ad=min_r C_r; C_max=max_r C_r; C_fix=capacity of the vertexwise minimum field; and C_upper=capacity of the vertexwise maximum field. Disconnected base-graph pairs have capacity minus infinity. C_fix <= C_ad <= C_max <= C_upper and C_ad <= E.

Only common thresholds in [0,E] have endpoints accepted in every scenario. Preserve the full exact threshold classification, not a selected grid:

| Label | Interval, intersected with [0,E] | Meaning |
| --- | --- | --- |
| fixed_route | (-infinity,C_fix] | One fixed route is accepted in every scenario |
| adaptive_routes_only | (C_fix,C_ad] | Each scenario admits a route, but no single route works in every scenario |
| mixed_scenario_connectivity | (C_ad,C_max] | Some scenarios connect and some do not |
| joint_obstruction_box_connected | (C_max,C_upper] | Every joint scenario is disconnected, while the vertexwise upper box admits a route |
| upper_box_obstruction | (C_upper,infinity) | Even the upper box has no route |

All five interval records are retained, including empty intervals. At exact equalities, use the closed upper boundary of the preceding category. The first interval includes zero when the base graph is connected. The all-scenario obstruction interval is (C_max,E], if nonempty; the box obstruction interval is (C_upper,E], if nonempty. These are threshold ranges, not confidence intervals or probabilities.

Store every endpoint's E, all 27 C_r values, four capacity summaries, all five intervals, and the two union-level obstruction intervals. Thus there are 69,012 scenario-pair capacity entries. Report category nonemptiness counts over all 2,556 pairs, with explicit denominators. A pair may have nonempty intervals in several categories, so those pair counts are not a partition. Preserve all native interval endpoints and widths. No aggregate width is interpreted as prevalence over thresholds. Primary statements concern common absolute thresholds only. Scenario-specific thresholds based on min(f_r(s),f_r(t)) are not evaluated in this branch and may not replace E.

## Unknown values and eligibility

A missing/undefined source scenario entry remains `None` and denotes an unresolved nonnegative score in [0,+infinity]. Do not replace it by zero as if observed. For computation only, construct lower and upper completion fields with zero and infinity respectively, preserving a flag for every unknown cell. Each monotone metric E, C_fix, C_ad, C_max and C_upper receives sound lower and upper bounds by evaluating these two completions.

No exact five-way classification is emitted for an incomplete pair-field collection. A threshold-specific helper may certify a label only if all compatible completions have that label, using lower/upper metric bounds and guaranteed endpoint acceptance; otherwise it returns unresolved. The helper also distinguishes certainly ineligible endpoints from unresolved endpoint eligibility. These are conservative sufficient certificates, not an exhaustive correlated-completion solver. Retain all 2,556 pairs with exact or bounded status. No unknown-dependent exclusion or chosen substitute scenario is allowed.

## Evidence and failure rules

The source ratios, fields and scenarios are already exposed. New summaries remain uncomputed until root approves the protocol and implementation hash. The mathematical derivation is elementary minimax quantifier accounting and bottleneck monotonicity; no novelty claim is made. The biological scope remains one selected minigene assay domain with dependent preparations.

Synthetic checks must verify known-answer parallel-route, series-route, mixed-scenario, endpoint-ineligible, tied, zero and disconnected examples. Compare the implementation against independently enumerated simple paths on every three-vertex undirected graph with every ordered pair of ternary fields. Exercise unknown values against every binary completion on a fixed three-vertex path, preserving exact boundaries. Test vertex/scenario permutations and a full 72-vertex, 27-field constant-score fixture to authenticate pair/scenario coverage. No biological source or result is opened by the test program.

Stop on an algebraic counterexample, source or schema mismatch, lost pair/scenario, false certificate, inconsistent tie rule, nonfinite defined score or unexplained numerical discrepancy. Finite verification supports the implementation within the stated range and does not replace the proof. Expansion, scenario-specific threshold analysis, correlations and new fitted models require a separate decision.
