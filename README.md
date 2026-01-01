# Synonymous RNA accessibility

Exact and thermodynamic analyses of synonymous sequence spaces under RNA pairing constraints.

Author: Vivaan Patni, Independent Researcher. Contact: contact.vpatni@gmail.com.

This repository contains the computational workflow, input provenance, exact verification cases, derived data, and reproducibility documentation. Manuscripts and journal submission materials are maintained separately.

## Scientific scope

Compatibility, folding-model admissibility, measured RNA structure, and organismal fitness are distinct. Disconnection is reported only for completely enumerated domains or a proved structural invariant. All counts use exact integers. Native mappings must match the reference assembly used in the source study. A sampled graph is never treated as proof of disconnection.

## Reproducibility

The analysis specification is in `analysis/specification.json`. Software versions, source checksums, exclusions, and validation results accompany derived data as the corresponding analyses finish. Install the pinned environment and run the tests before reproducing results. Third-party data are retrieved from their original or version-pinned sources rather than silently relicensed.

The initial specification is historical. The actual primary five-codon eligibility rules are frozen in `analysis/native_specification.json`, with the measured-coverage correction in `analysis/eligibility_correction.json`. Expansion controls, the discovery boundary and the prospective validation procedure have separate versioned specifications in `analysis/expansion/`. Previously examined outcomes are not relabeled as prospective validation.

### Environment

The tested runtime is CPython 3.12.14. Install `requirements-tested.txt` in an isolated environment and run:

```text
python -m pytest tests -q
python analysis/exact/run_final_audit.py --out /path/to/fresh-register-check
python analysis/exact/check_automata_certificates.py
python analysis/expansion/verification/run_verification.py --outdir results/expansion/fresh-general-check
python analysis/expansion/verification/verify_feasibility_2sat.py --outdir results/expansion/fresh-feasibility-check
```

The general-graph and feasibility drivers require fresh output directories. Assertions must remain enabled; optimized Python executions are rejected by the verification drivers. Source hashes are recorded before and after the general/feasibility checks. Finite verification supports the stated proofs and never proves an unstated universal claim. Numerical root estimates supplement exact integer-vector spectral certificates.

R 4.6.1 was used for the primary 10,000-resample transcript-cluster bootstrap. Its script uses base R only. The C++23 reference kernel was compiled with GCC 16.1.0; it is optional for the Python verification suite. RNAstructure 6.6 is an external sensitivity dependency, distinct from ViennaRNA, and its executables and parameter files are not redistributed.

Recover the complete primary input maps and verify the frozen cohort with the pinned environment active:

```text
python analysis/prepare_native_inputs.py --output /path/to/input-cache --verify-frozen
```

Add `--cached-raw /path/to/weeks194_orig` to use existing checksum-verified raw files without network access. The script verifies the bundled U00096.2 GenBank, 194 exact unique transcript mappings, 409 CDS occurrences, and all-strand coding overlap exclusions. It reconstructs the full 393,867-row coding-position map. Verification runs selection in a temporary directory and requires exact equality of all 404 frozen blocks and their complete exclusion audit, without overwriting frozen data. The completed check is recorded in `results/native_input_reproduction.json`.

### Native computations

The archive includes the full compatible genotype sets, point-mutation edges, retention fields, filtered components and mapped positions. Recomputing folding is more expensive than replaying graph summaries:

```text
python analysis/native_models.py run --workers 8
python analysis/summarize_native.py
Rscript analysis/native_inference.R .
Rscript analysis/native_inference.R . deduplicate
python analysis/expansion/conditional_vertex_null.py --workers 8
python analysis/expansion/retention_filtration.py --workers 8
python analysis/expansion/smoothness_controls.py --workers 8
python analysis/expansion/additive_score_control.py --workers 8
python analysis/expansion/context_filtration.py
python analysis/expansion/extension_inference.py
python analysis/expansion/residual_stress.py
```

The original folding runner reuses existing native result files. To recompute it, direct a new checkout/result copy to a fresh result location or move the existing native outputs aside first; do not overwrite the only preserved receipt. The new score-control runners use content fingerprints for method source, inputs, parameters and runtime versions, and invalidate stale caches. The bootstrap is descriptive uncertainty across the modeled transcript collection, not biological-replicate or evolutionary uncertainty.

Independent original-PAIR-MaP source provenance and reproduction instructions are in `analysis/pairmap/README.md`. That source was already examined during discovery. The first 24 primary windows evaluated with RNAstructure, context variants and larger complete blocks are sensitivity analyses, not new independent biological validation. Larger-window gates and every resource failure are retained in `results/expansion/longer_selection.json`.

### Interpreting the score controls

Uniform score permutations retain component-specific score multisets and destroy neighbor correlation. Degree conditioning additionally retains coarse graph-position composition at cutoff 0.5. The exact small-space control conditions on an edge-roughness interval. Heat-rank fields preserve score multisets and approximately match mean local roughness; they are a generic reference model, not uniform conditional permutations. Additive fits are descriptive approximations using the same observed score fields. Their in-sample fit is not predictive validation.

Generic smoothness attenuates approximately 94–95% of the shuffled-score contrast. Additive approximations produce slightly less fragmentation than the folding fields. Consequently the archive supports exact compatibility geometry and model-derived score barriers, not an RNA-specific evolutionary connectivity-preservation law. `analysis/expansion/claim_evidence_map.md` separates simultaneous compatibility, alternative conformations and physical ensembles.

### Frozen computational prediction pilot

`analysis/expansion/barrier_pilot_freeze.json` and the complete predictions in `results/expansion/barrier_pilot/predictions.json` were committed as `01c3243` before full-ensemble evaluation. This is a separate synthetic computational protocol, not a replacement for independent biological validation. All 24 scaffolds, 72 landscapes, 48 interventions, matched controls and failures remain available. The isolated-hairpin approximation misses five threshold-defined barriers, has worse lifetime error than zero prediction, and succeeds on only 12/24 covered intervention directions. Every three-state additive fit is saturated. The pilot fails its expansion gates; no functional rescue or folding-specific engineering claim follows.

```text
python analysis/verify_prediction_pilot.py
python analysis/expansion/barrier_prediction_pilot.py validate
python analysis/verify_bundle.py --out /path/outside/repository/integrity.json
```

The prediction table is `results/expansion/barrier_pilot/prediction_table.csv`. It is a model-test table, not an experimentally validated design list. The `prepare` stage refuses to overwrite frozen predictions; do not delete that freeze to obtain a more favorable candidate set. `validation_source_audit.md` records the specific missing structure/probing export that prevents closure of the independent biological-material branch. Discovery-window extensions and independent engine implementations remain sensitivity analyses.

### Exact uncertainty and additive product queries

The current programme separates full marginal boxes, joint affine uncertainty and fixed versus scenario-dependent routes. Sharp interval bounds are deterministic identified sets, not confidence intervals. `analysis/programme/product_closure_proof.md` proves a direct additive-product endpoint query and its witness bound. Additive merge-tree composition is established prior work; `analysis/programme/product_prior_manifest.json` records the unresolved priority comparison with the closest detailed construction. No first or fastest algorithm claim is made.

```text
python analysis/programme/check_identification.py --out work/fresh_identification.json
python analysis/programme/joint_uncertainty_verify.py --outdir work/fresh_joint_uncertainty
python analysis/programme/product_closure_checks.py --output work/fresh_product/results
python analysis/programme/product_closure_independent.py --outdir work/fresh_product_independent
python analysis/programme/forest_depth_regression.py
python -m pytest tests/test_product_closure.py tests/test_forest_extensions.py tests/test_overlap_quotient.py -q
```

The product checks require a new output location. The independent audit uses explicit product traversal and a separate factor closure. Its present result is in `results/programme/product_closure_independent_release`. The forest regression compares a 1,200-event height-two path against the independent Entringer recurrence. The authoritative overlapping-obligation API receipt is `results/programme/overlap_v2_verified_api`; earlier source hashes are historical, not verification of the amended API.

### Acquisition validation remains a failed efficiency result

The frozen held-out acquisition benchmark contains 3,800 worlds and 22,800 policy runs, with every result retained. The proposed path/cut rule fails the predeclared efficiency hurdle and uses 0.398 more queries on average than the strongest forward baseline. The 95% seed-block interval [0.261, 0.535] describes Monte Carlo variation in this fixed synthetic design. The examined validation set is retired from tuning. See `results/programme/acquisition_v1/validation_report.md` and its frozen protocol for the exact paired design, all comparators and failures.

### Measured RNA pilot

A complete 32-genotype reference-protein fibre from the Gill2026 deposit was frozen before reading RNA outcomes. All predictions were committed before the held-out RNA library. The pairwise predictor fails its primary hurdle: RMSE 0.551 versus 0.466 for the constant baseline, with failure at both pseudocount sensitivities. This is conditional reproduction of a held-out RNA library on the same genotypes. Biological independence, RNA/DNA pairing and calibrated measurement uncertainty remain unresolved. No biological interaction or folding mechanism is inferred from positive observed-field barriers. The complete protocol, source attribution, chronology and replay instructions are in `results/programme/gill_reference_pilot/README.md`.

### Figure-source reproduction

Panel data are in `results/figure_data/`. The fresh rendering source writes assets only outside the computational repository:

```text
python analysis/render_figure_suite.py --out /path/outside/this/repository/figures --figures 1,2,3,4,5,6,S1,S2,S3,S4
python analysis/render_identification_figure.py --out /path/outside/this/repository/figures
python analysis/render_joint_uncertainty_figure.py --out /path/outside/this/repository/figures --qa work/figure8qa
python analysis/render_gill_pilot_figure.py --outdir /path/outside/this/repository/figures --qa /path/outside/this/repository/figure9qa
```

Figure 9 uses the frozen panel-data manifest `results/figure_data/Figure_9.json`. Its renderer verifies the saved training and validation source hashes, recomputes the displayed summaries, and exports the complete failed prediction comparison without opening the source workbook or refitting models. The required `--outdir` and `--qa` arguments are confirmed by `python analysis/render_gill_pilot_figure.py --help`.

No manuscript, submission files, article PDFs, figure exports or publishing automation belong in this archive. The continuous-integration workflow only verifies computations.

### Commit privacy

Local commits were made through git-privacy 2.3.0 with `privacy.pattern=Mdhms`, UTC dates and author/committer dates reduced to the start of their year. `requirements-privacy.txt` pins the tested privacy tooling separately. Commit-date sanitization does not conceal GitHub push events or public release timestamps.
