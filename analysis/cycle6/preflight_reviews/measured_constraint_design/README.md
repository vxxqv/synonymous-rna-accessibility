# Bounded mathematical decision handoff

- `report.md`: source-level priority decision, sparse separator and necessary-detour contract, conservative marginal uncertainty, counterexamples and FAS gate-aware applicability.
- `prior_map.md`: targeted original-source overlap map and explicit extraction/access limits. No novelty claim.
- `ordinal_interval_characterization.md`: exact interval/order equivalence, including finite lower/+infinite upper bounds, ties, constants and rational cancellation certificates.
- `protocol_math_review.md`: issues and inspected correction of the root-owned mutable draft. It does not approve the changing production implementation.
- `fixtures_v2.json`: authoritative finite structural examples; `fixtures.json` is the earlier cube-only receipt.
- `ordinal_interval_checks_v3.json`: authoritative 2916-instance interval reference receipt. The earlier finite-only and v2 receipts are historical.
- `receipt.json`: final file/input hashes and scope.

## Reference API

`ordinal_interval_audit.py` uses the Python standard library and exact `Fraction` arithmetic. Import is safe and writes no output. Functions are:

```text
inequalities(X, intervals) -> (forced_pairs, rows)
solve(rows, dimension) -> rational coefficient tuple or None
fitted_values(X, intervals, coefficients) -> (latent_scores, assigned_responses)
```

Each row is `(coefficient_tuple, right_hand_side)` for `a*x>=rhs`. Upper bound `None` means +infinity; lower bounds must be finite. Inputs are mathematical rational values. This is a small reference and synthetic oracle, not a scalable production LP solver; Fourier-Motzkin elimination can grow exponentially. For a large observed result, independently checking an exact supplied primal/dual certificate does not require rerunning this elimination.

Copy the script alone to a fresh working directory and run `python -B ordinal_interval_audit.py` to create a new `ordinal_interval_checks_v3.json`. Existing output and optimized Python are refused. `fixtures.py` similarly creates `fixtures_v2.json` and refuses overwrite. No Git, network, assays, optional packages, author functions or raw sources are used by these scripts.

All outcomes remain unread by this lane. The metadata source report did contain the already-disclosed five illustrative source-README occurrence rows, which remain part of the source exposure ledger; they were not used here. None is a fresh protected phenotype or a selected validation result. No shared file, repository, manuscript, commit or external service was changed.
