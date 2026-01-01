# Synonymous RNA accessibility

Exact and thermodynamic analyses of synonymous sequence spaces under RNA pairing constraints.

Author: Vivaan Patni, Independent Researcher. Contact: contact.vpatni@gmail.com.

This repository contains the computational workflow, input provenance, exact verification cases, derived data, and reproducibility documentation. Manuscripts and journal submission materials are maintained separately.

### Current measured-function feasibility result

`analysis/cycle5/README.md` documents the frozen joint-function pilot on 47 intended-synonymous four-genotype designs. Its protocol, explicit uncertainty amendment, global target masks, independent implementation checks and source-provenance limits are committed before new outcomes are accessed. The source archive remains external to Git and is checksum pinned. This is an author-preprocessed within-library holdout, not independent biological replication. All 12 deposited-series scenarios, source-mean gates, failed prediction criteria, observed-component bounds and group influence diagnostics are prescribed in advance. Original draft states and the initially incomplete diagnostic revision remain as historical records; `cambray_pilot/protocol_v1.md` and `freeze_v1.json` identify the authoritative specification.

`analysis/cycle5/CURRENT_RESULTS.md` gives current portable verification commands. The completed pilot is recorded in `results/cycle5/pilot_closure_report.md` and `pilot_decision.json`. Reserved coverage is six squares in six groups, below the frozen ten/eight minimum. One reserved finite square passes all twelve series scenarios and a fixed-low-end empirical-range certificate, but biological uncertainty is uncalibrated. Both continuous prediction screens fail. Both convex mean classification baselines are structurally negative in exact square geometry. Every failed comparison, missing target and group deletion is retained. Both selected arms are now exposed and unavailable for new prospective tuning.

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

Generic smoothness attenuates approximately 94 to 95% of the shuffled-score contrast. Additive approximations produce slightly less fragmentation than the folding fields. Consequently the archive supports exact compatibility geometry and model-derived score barriers, not an RNA-specific evolutionary connectivity-preservation law. `analysis/expansion/claim_evidence_map.md` separates simultaneous compatibility, alternative conformations and physical ensembles.

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

### Multidimensional correlated uncertainty

`analysis/programme/affine_polytope/` extends the deterministic uncertainty analysis to bounded rational parameter polytopes. Its reviewed proofs identify exact ranges at fixed latent dimension, establish variable-dimension complexity boundaries, and characterize when one fixed route is as robust as scenario-dependent routes. Generic lattice-polynomial and arrangement methods are established prior; the module makes no first-algorithm claim and supplies no calibrated biological confidence region.

```text
python analysis/programme/affine_polytope/verify_exact.py --out work/fresh-affine-polytope
```

The portable exact replay compares 82 cases against a separately authored path and separator oracle, checks 400 rational quantities, and retains the original corpus and disclosed fixture correction. The README in that directory gives the precise proof, review and provenance boundaries. These additions postdate the original version 1.0 archive and must not be attributed to that archive until a later verified release includes them.

### Bounded interaction certificates

`analysis/programme/interaction_certificates/` contains a portable, byte-preserving replay of the reviewed projection and path bounds. Small blocks can certify an exact interacting barrier in some cases, but even one interaction edge can defeat singleton projections. A switch reduction establishes an exact-capacity hardness boundary for bounded sparse interactions, conditional on the cited established MaxCut theorem. No novel elimination method or general efficient interacting-path algorithm is claimed.

```text
python -B analysis/programme/interaction_certificates/replay.py --out work/fresh-interaction-certificates
```

The module preserves the original proof, independent review, corrected cost contract and accounting erratum. Read `revision_2/claim_and_proof.md` and `revision_2/complexity_contract.md` in the extracted payload as the current statements. Witness storage, retained envelopes, actual elimination order and baseline-component restriction are explicit. The replay needs no Git or network and does not reauthenticate historical Git events. These additions postdate the published version 1.0 archive.

### Acquisition validation remains a failed efficiency result

The frozen held-out acquisition benchmark contains 3,800 worlds and 22,800 policy runs, with every result retained. The proposed path/cut rule fails the predeclared efficiency hurdle and uses 0.398 more queries on average than the strongest forward baseline. The 95% seed-block interval [0.261, 0.535] describes Monte Carlo variation in this fixed synthetic design. The examined validation set is retired from tuning. See `results/programme/acquisition_v1/validation_report.md` and its frozen protocol for the exact paired design, all comparators and failures.

### Exact observation costs and optional batches

The bounded audit in `analysis/programme/observation_cost/` uses standard finite decision analysis on six hidden vertex states. All 384 cost/prior/action-set cases and all complete decision tables are retained. Optional batches improve the exact optimum in 184 of 192 matched comparisons, while the cheapest-per-unseen-bit rule is never optimal on this grid. These exact conditional expectations are not measured assay costs, general policy guarantees or a repair of the failed held-out acquisition result.

```text
python -I -B analysis/programme/observation_cost/replay.py --out work/fresh-observation-cost
```

The replay requires the pinned preinstalled NetworkX 3.7, uses no runtime network or Git, and compares all scientific outputs to their independently checked receipts. The original proof is accompanied by a clarification for zero-probability histories under generic correlated priors. The frozen experiment uses full-support product priors and its results are unchanged. This addition postdates published version 1.0.

### Measured RNA pilot

A complete 32-genotype reference-protein fibre from the Gill2026 deposit was frozen before reading RNA outcomes. All predictions were committed before the held-out RNA library. The pairwise predictor fails its primary hurdle: RMSE 0.551 versus 0.466 for the constant baseline, with failure at both pseudocount sensitivities. This is conditional reproduction of a held-out RNA library on the same genotypes. Biological independence, RNA/DNA pairing and calibrated measurement uncertainty remain unresolved. No biological interaction or folding mechanism is inferred from positive observed-field barriers. The complete protocol, source attribution, chronology and replay instructions are in `results/programme/gill_reference_pilot/README.md`.

### Additional protein-preserving measurement groups

The multiclass Gill test contains 1,124 development groups and 1,052 validation groups from one reporter. The previously examined reference-protein stratum is quarantined. The frozen normalized residual cosine is 0.07788 against the declared 0.5 target; residual prediction increases mean squared error by 47.7% over zero. High RNA-only and DNA-only agreement does not establish normalized residual reproducibility. All pseudocounts, 45 normalization-pair comparisons and single-stratum omissions are retained. Groups share one reporter and unresolved preparation identities. No independent-gene, biological-noise or folding-mechanism inference follows.

```text
python analysis/programme/retrospective/gill_multiclass_v1/independent/replay.py --repo-root . --output /new/path/outside/repository/gill-replay
```

The offline independent replay verifies 3,352,385 comparisons from archived processed counts. Its guards reject optimized Python, existing outputs and outputs inside the repository. Historical access order is documented separately and is not authenticated by numerical replay.

### Full-length reporter source-family transfer

The selected PERSIST cohort has 22 development and 15 held-out constructs encoding one protein in one context. Protocol and implementation were frozen before features and selected outcomes. All six prediction vectors were committed before validation. The primary model adds one explicitly recalculated mean unpaired probability (MUP) to composition and codon features. It fails the 10% improvement target: MAE 0.031442, compared with composition alone at 0.030139. The matched codon comparator improves by only 0.119% on adding MUP. Every held-out row and every omission diagnostic is retained; no model or subset was retuned.

```text
python -I -B analysis/programme/retrospective/persist_v2/independent_replay/replay.py --repo . --out /new/path/outside/repository/persist-replay
```

This independently reconstructs tuning, predictions, metrics and influence diagnostics from immutable features and the 37 selected target cells. It treats MUP as an input and does not refold sequences. Deposited cellular coefficient units are not calibrated physical half-lives; source-family separation does not establish independent biological preparations or complete mutation paths. The protocol's original Windows binary hashes remain historical provenance. See the module README for scope and source attribution.

The separate frozen feature replay reconstructs all 83 features for all 37 selected constructs, including the fixed ViennaRNA ensemble calculation. Full Windows and Linux replays reproduce every value exactly. The initial Linux guard failure was traced to CRLF versus LF in an otherwise identical parameter file; the original failure and a narrowly scoped line-ending correction remain archived. It uses the original scientific functions and parameters, checks selected sequence hashes, and stops on warnings or discrepancies. The archived original wrapper and its pre-execution warning-guard correction remain visible. This verifies the specified calculation, not biological validity or measured structure.

```text
python -I -B analysis/programme/retrospective/persist_v2/feature_replay/check_guards.py --output /new/path/outside/repository/persist-feature-guards
python -I -B analysis/programme/retrospective/persist_v2/feature_replay/replay.py --repo-root . --output /new/path/outside/repository/persist-features
```

These commands require the pinned dependencies already installed. They use no Git, network, outcome access, refitting, retries or fallback settings. The runtime records its actual platform and binary hashes; original Windows binary identity is not required for cross-platform replay.

### CNOT3 association and baseline lifetime

The original published-cohort gate remains failed: independent readers recover 5,480 rather than 5,481 genes, and the historical missing identifier remains unresolved. One documented operational-cohort amendment preserved all scientific model choices. It is exploratory after source-outcome exposure, not untouched or independent experimental validation. No original value was corrected to repair the gate.

The fixed models use 4,304 training and 1,112 held-out genes with related-protein grouping. Primary MAE is 0.339316. Association improves MAE by only 0.113113% relative to codon and 0.041944% relative to the strongest primary baseline, failing the 5% criterion. All nine models, fixed sensitivities, component deletions and 10,008 held-out prediction/error rows remain. Independent numerical reconstruction passed 59,268 checks over 240,144 numerical values. The examined holdout is retired from tuning.

`results/programme/cnot3_v2_exploratory/result_report.md` gives complete interpretation and source-exposure limits. The original failed gate and source evidence are preserved alongside the explicitly amended successor. Half-life uncertainty is uncalibrated; component deletions are influence diagnostics. Elementary kinetic bounds are conditional identification sets, not confidence intervals or a pathway-specific molecular rate estimate.

```text
python -I -B analysis/programme/retrospective/cnot3_v2_exploratory/source_diagnosis/portable_replay/replay.py --repo . --out /new/path/outside/repository/cnot3-source
```

This byte-pinned offline source replay checks the complete deposited population and arithmetic diagnosis without Git or network. The final independent numerical replay and its ten guard tests are available in `analysis/programme/retrospective/cnot3_v2_exploratory/final_numerical_replay/`. Run the following from the archive root with the pinned runtime:

```text
python -I -B analysis/programme/retrospective/cnot3_v2_exploratory/final_numerical_replay/replay.py --archive-root . --out /new/path/outside/repository/cnot3-numerical
```

It verifies every frozen prediction and all scored diagnostics without importing author code. The source and numerical replays require new output directories. Historical access order is not authenticated by numerical reconstruction.

### Figure-source reproduction

Panel data are in `results/figure_data/`. The fresh rendering source writes assets only outside the computational repository:

```text
python analysis/render_figure_suite.py --out /path/outside/this/repository/figures --figures 1,2,3,4,5,6,S1,S2,S3,S4
python analysis/render_identification_figure.py --out /path/outside/this/repository/figures
python analysis/render_joint_uncertainty_figure.py --out /path/outside/this/repository/figures --qa work/figure8qa
python analysis/render_gill_pilot_figure.py --outdir /path/outside/this/repository/figures --qa /path/outside/this/repository/figure9qa
python analysis/render_gill_multiclass_figure.py --outdir /path/outside/this/repository/figures --qa /path/outside/this/repository/figure10qa
python analysis/render_persist_figure.py --outdir /path/outside/this/repository/figures --qa /path/outside/this/repository/figure11qa --manifest /path/outside/this/repository/Figure_11.json
python analysis/render_cnot3_figure.py --data results/figure_data/Figure_12.json --expected-data-sha256 3f40c1d627effd1fb0934b6ff676952bcca7030ae88e00b774b5576d061dbf31 --out /new/path/outside/repository/figure12
```

Figure 9 uses the frozen panel-data manifest `results/figure_data/Figure_9.json`. Its renderer verifies the saved training and validation source hashes, recomputes the displayed summaries, and exports the complete failed prediction comparison without opening the source workbook or refitting models. The required `--outdir` and `--qa` arguments are confirmed by `python analysis/render_gill_pilot_figure.py --help`.

No manuscript, submission files, article PDFs, figure exports or publishing automation belong in this archive. The continuous-integration workflow only verifies computations.

### Commit privacy

Local commits were made through git-privacy 2.3.0 with `privacy.pattern=Mdhms`, UTC dates and author/committer dates reduced to the start of their year. `requirements-privacy.txt` pins the tested privacy tooling separately. Commit-date sanitization does not conceal GitHub push events or public release timestamps.
