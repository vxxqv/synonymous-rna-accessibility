# Reproducing the measured-data feasibility gates

This cycle did not complete a biological prediction model. Both URA3 executions intentionally stopped when their source-data guards failed. The mitochondrial ensemble comparison stopped at its pre-outcome coverage gate. Read `measured_function_closure_v6.json` before interpreting any candidate mapping as usable measurements or biological evidence.

Use the project's tested Python environment (Python 3.12.14, NumPy 2.5.4, Biopython 1.88). The fetcher uses only Python's standard library. No author scripts, R packages, folding engine or raw-read pipeline is needed for these gates.

From the repository root, choose an ignored source directory and fresh output locations:

```text
python analysis/programme/fetch_measured_sources.py --work data/raw/measured_function
python analysis/programme/measured_data_audit.py --work data/raw/measured_function --out work/measured_audit_replay
python analysis/programme/verify_measured_data.py --work data/raw/measured_function --out work/measured_verification.json
```

`data/programme/measured_source_manifest.json` identifies eight core external inputs by public URL and SHA-256. `measured_endpoint_source_manifest.json` adds the original single-mutant supplement and two GEO metadata files for direct source reconciliation and sample-ID checks. External author code and papers are not bundled. The only accepted text preparation is explicit terminal-newline/CRLF normalization for FASTA files, and even then the resulting file must match the recorded hash. Other source changes fail verification. Source availability or a hash mismatch must be reported, not silently replaced with a newer dataset. The table containing mitochondrial functional measurements is fetched to verify provenance/reference compatibility, but the analysis does not read its numerical phenotype values.

The genotype audit checks only text labels and numeric-cell presence. Its 6470 condition groups are **candidates**, not 6470 independent or usable measurements. Zero-coded fitness cells are not interpreted as missing, neutral or lethal. All four CDS translations and sequence hashes are independently replayed using a literal genetic-code table. `independent_verification.json` records 3449 distinct reconstructed coding genotypes and eight source-hash checks.

The structural audit checks all eight population-profile sequence strings against the native reference with its N3107 removed. The function reference is checked against the unmodified native reference. Local cluster sequences are anchored to genes; plotting coordinates are not used as genomic coordinates. The frozen 30-nt cluster-feature gate has six of eleven units, below the required nine. The 15/45-nt windows have six/seven units. Do not fit the six-unit subset or impute cluster rates.

The original pilot stopped at a zero-valued single. A separately named correction honors the originally declared positive-finite eligibility, and includes a separately frozen guard for values below the described 0.5 source cutoff. To reproduce the second stop:

```text
python analysis/programme/ura3_function_pilot_v2.py --work data/raw/measured_function --audit results/programme/measured_data_audit_v1 --out work/ura3_guard_replay
```

It must exit nonzero at J11262, G321T/C336T, reported fitness 0.326463358865135. This is the recorded source-consistency outcome. It is **not** a biological hypothesis rejection or evidence of an accessibility barrier. No prediction seal, fitted model or validation metrics should be produced. Do not remove the guard to manufacture a model result. The existing stopped-run manifests retain commands, source hashes, synthetic controls and errors; absence of declared prediction results is expected and explicitly recorded.

The single-mutant source-reconciliation audit is described in `endpoint_reconciliation_v4.json` and `endpoint_decision_v5.json`. It compares the current report with the original S4 and the already-published C315T example; it does not read the planned validation doubles. Exact zeros and the subcutoff row cannot be individually assigned a historical exclusion cause from the released processed files. Recover genotype-by-replicate counts/fitness, barcode maps, WT denominators, growth generations and exclusion flags/order before reconsidering the gate.

All decisions are chronological records. `measured_function_protocol_v1.json` fixes the assay, regional split, endpoints, comparisons and uncertainty scope. Later files state source problems and implementation amendments without rewriting that freeze. File hashes remain unchanged by the separately requested author-attribution rewrite; `data/programme/vxxqv_attribution_map.json` maps historical commit identifiers to their content-identical author-corrected commits. Dates remain coarsened through git-privacy.
