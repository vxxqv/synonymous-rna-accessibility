# Full operational candidate-universe check

The remaining independence gap is closed. A separate check independently reads only A/I/J/K across all22,848 data rows, before comparing membership with the frozen selector. It finds exactly7,359 finite-IP candidates and15,489 exclusions. Every excluded row lacks a finite value in all three I/J/K fields; no Decimal-versus-float finiteness disagreement or ID collision was found. The7,359 row set exactly matches the original exposure selector, and every source GeneID matches the authenticated earlier GeneID-only extraction.

Together with the preserved independent B/C/D audit on those candidates, this authenticates the full **operational5,480-gene cohort** under the frozen rules. It rules out a candidate missed by the original stage1 reader within the source table. It still does not recover the publication's5,481-member universe or establish a published typo.

The check uses the same previously independently authored, hash-pinned SAX collector, configured in its new process to select only A/I/J/K. Original helper bytes and `observed_v2` are unchanged. No B/C/D outside the frozen candidates, F/G, other outcomes, grouping or models were opened or run. These A/I/J/K cells had already been exposed by the original authorized reproduction; this narrow authorization adds no new endpoint exposure.

From the workspace root, reproduce into a fresh checkout of this review folder before its receipt exists:

```text
python -I -B work/decay_cycle4/source_discrepancy_review/check_candidate_universe.py
```

The fixed output is `candidate_universe_check.json`; overwrite is refused. It records all input/helper hashes before and after, exact candidate-set hashes and the scope. Receipt SHA256: `8e2788e581f5729af27b1ca0bcc473b3233792dee82a9608f75a94bfc49667b5`.

The preceding bounded half-life arithmetic can be replayed independently into a fresh owned directory:

```text
python -I -B work/decay_cycle4/source_discrepancy_review/audit_source.py --out work/decay_cycle4/source_discrepancy_review/fresh_replay
```

Neither replay needs network, Git, installation, an original absolute Windows path or any new outcome workbook. They require the unchanged source workbook and existing whitelisted inventory/exposure/frozen-result files in the recorded workspace-relative layout. The earlier source manifest and reports remain historical evidence; this addendum closes only their explicitly identified candidate-universe limitation. The failed gate remains unchanged.
