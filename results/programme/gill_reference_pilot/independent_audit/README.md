# Independent Gill result audit

`portable/audit.json` and `portable/manifest.json` are authoritative for the current auditor. The Mathbox schema-2 manifest validates with repository root set to this computational archive. `current/` records the preceding implementation and is historical.

The initial `audit.json`, plus the source/result copies under `history/`, retain the first finite result checks. The current auditor adds an explicitly requested source-only replay mode for archives without Git history. Its scientific arrays, counts and differences were verified identical to the earlier audit and the default full-history mode. No model, count, endpoint or result was changed.

`preoutcome_implementation_archive.json` preserves the historical synthetic-only implementation audit and records its archived source path changes. That archive was not replayed on observed outcomes.

To reproduce from the repository root, use the audited Python environment and run:

```
python analysis/programme/retrospective/gill_result_independent.py --outdir work/gill-independent-reproduction
```

Choose a fresh output directory; existing results are never overwritten. The default requires a full Git clone with the frozen protocol and prediction commits. For a source archive, add `--source-only`; this retains all arithmetic and stored-hash checks but explicitly reports that historical Git evidence was not verified. The auditor opens saved count/prediction JSON files and, in default mode, Git blobs. It never opens the observed workbook or raw reads. Complete instructions and limitations are in the parent README.
