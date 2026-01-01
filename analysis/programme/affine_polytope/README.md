# Correlated uncertainty on finite accessibility graphs

This module identifies deterministic capacity and barrier ranges when vertex scores depend affinely on shared parameters in a bounded rational polytope. It distinguishes a route available in some scenario, a route available in every scenario, and one fixed route available in every scenario. It supplies no confidence level for a biological measurement model.

The reviewed statements are in `claims.md`:

* Exact score-order enumeration is polynomial for each fixed parameter dimension, with the exponent depending on that dimension. All preprocessing is included in the stated arithmetic bound. Rational bit costs are additional.
* When parameter dimension varies, possible threshold connectivity is NP complete and guaranteed threshold connectivity is coNP complete under the explicit constructions and representation given there.
* Fixed and scenario-dependent guaranteed capacities agree for every compact score uncertainty set exactly when the graph has one inclusion-minimal connecting vertex set. Uniqueness of a simple path is a stronger condition than necessary.

Weighted lattice-polynomial normal forms and score-order representations are established prior. The first statement applies that machinery to accessibility; the third is a structural consequence. Precise priority of the complexity encodings remains unresolved. No first or fastest algorithm claim is made.

Run from the repository root with Python 3.12:

```text
python analysis/programme/affine_polytope/verify_exact.py --out work/fresh-affine-polytope
```

The destination must not exist. The check uses only the Python standard library. It compares the arrangement implementation against a separately authored exact path and separator oracle using Fourier-Motzkin elimination. It also checks all 1,098 labelled simple graphs on two through five vertices and 259 specified Boolean formulas. Finite agreement supports the supplied proofs within the checked ranges and is not itself a proof.

The active corpus contains 82 cases and 400 exact rational comparisons. All 81 original cases and their values are preserved. One original fixture was misleadingly named: its capacity is `abs(theta1+theta2)`, and its extrema already occur at parameter corners. The corrected name is `rank_one_sum_two_parameters`. A separately disclosed six-vertex example has capacity `max(abs(theta1),abs(theta2))`; every square corner misses its minimum at the origin. The original inputs, name and receipt remain in `historical_v1/`. The change ledger and independently authored review document this correction.

The review was separately authored with inherited project context, not blinded. `reference/` contains the unchanged review code, argument and receipts. Original receipts retain historical workspace-relative paths and their original source hashes. The command above is the portable entry point. `integration_sources.json` records the byte-for-byte source copies; `receipts/` preserves the corrected replay. No manuscript or figure exports are included.
