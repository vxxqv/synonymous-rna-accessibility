# Exploratory successor A: exact ordinal exclusion robustness

This is an adaptive follow-up after all Cycle 6 outcomes were exposed. It is not fresh validation. The primary scientific freeze, all primary results and the zero primary barrier count remain unchanged. Freeze this contract and tested implementation before computing the new robustness quantities.

## Question and input

How far must the complete empirical envelope be expanded before some codon-additive latent score and some arbitrary common nondecreasing response fit all declared observations? Use exactly the 72 Cycle 6 vertices, five codon factors, 27 scenarios and source result SHA-256 `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466`. No new source, genotype, threshold, response family or outcome-selected subset enters.

The native reference (index 46) equals one by definition and remains fixed at [1,1]. Perturbing that identity would not represent the chosen reference-relative measurement. Let a_v=0 for that reference and 1 otherwise.

## Two explicitly different metrics

1. Additive expansion epsilon >= 0: [max(0,L_v-a_v epsilon), U_v+a_v epsilon]. The result is in relative-enrichment units.
2. Symmetric multiplicative expansion gamma >= 1: [L_v/gamma^(a_v), U_v gamma^(a_v)]. The result is a factor relative to the reference. This is an existential envelope sensitivity, not an error confidence bound or automatically a simultaneous correction of all observed preparations.

A fixed persistent recovery-factor interpretation requires its own observation-model assumptions. It must not be confused with allowing a new correction in every preparation. No empirical calibration currently establishes a permissible bias bound.

## Exact calculation and decision

For additive expansion an original strict order i<j disappears at epsilon=(L_j-U_i)/(a_i+a_j), unless its denominator is zero. For multiplicative expansion use t=gamma^2. With exponent sum one the breakpoint is (L_j/U_i)^2; with sum two it is L_j/U_i. Positive orders above a zero upper score, or involving two fixed values, can be permanent. Ties remove all corresponding strict orders simultaneously.

The feasible set grows monotonically. Search the finite exact rational breakpoints for the first feasible system. Verify an exact dual contradiction whose supported orders remain active at every smaller parameter, and an explicit feasible latent vector and common response at the boundary. A floating solver status alone never counts. With nonnegative scores, squared multiplicative intervals [L_v^2/t^(a_v), U_v^2 t^(a_v)] are rational and preserve all strict orders. A rational feasible response in these squared units gives the original response by its nonnegative square root. Report an exact rational t and a rational enclosure of gamma. Treat permanent contradictions as an infinite radius only with a verified permanent-order dual.

Success is a sharp exact threshold or a certified interval/unbounded result, with both sides independently checked. Failure is an unresolved certificate or model/implementation ambiguity, not a favorable biological outcome. A small radius must be reported as fragility. No radius is called biologically robust without an independently supported error or recovery bound smaller than it.

## Complete influence diagnostics

- Omit each of the 72 observation rows in turn at the original envelope and report every exact feasibility result and certificate. Keep all other rows. Omitting the reference observation does not erase its role in the original measurement normalization. A genotype is indispensable only if deleting its observation makes the remaining system feasible.
- Omit each of nine output preparations and each of three input preparations from the scenario set, recompute the remaining coordinatewise envelope and calculate both complete-system radii. Report all twelve results. These sets are nested inside the original envelope, so deletion cannot remove the original contradiction and cannot reduce the radius. This monotonicity is an interpretation constraint, not evidence from independent validation.
- Preserve the selected five-comparison margins as historical witness diagnostics. They are not substituted for the optimized full-system radius.

No alternate normalization or source-reconstruction rule is introduced. The already computed 27 scenarios and corrected point remain the complete declared source-supported sensitivities. A different biological calibration would require new authenticated evidence and its own protocol.

## Verification and scope

Known-answer controls cover an accepted constant field, a reciprocal-sign square with exact additive radius 3/8 and squared multiplicative radius 4, a fixed-reference duplicate-row example, tied breakpoints, zero upper scores, permanent contradictions and an unknown interval. Independent review must verify the boundary proof, anchor treatment, all exact output certificates and the full influence inventory.

The generic radius routine handles unknown upper intervals. The observed-preparation adapter is explicitly restricted to the authenticated complete 72-by-27 primary input and stops if any entry becomes unknown; it never deletes that entry. Binary perturbation masks must be integers. Square-root enclosures use exact rational bounds and exact decimal strings, avoiding floating rounding claims. Pre-outcome review failures and their fixes are retained separately.

The method uses established order representation and finite active-constraint changes. No foundational novelty is claimed. No new molecular mechanism, biological selection, predictive generalization or barrier-specific interaction is inferred. Stop after this complete finite sensitivity calculation and independent review, then interpret it alongside the independently audited observation model.
