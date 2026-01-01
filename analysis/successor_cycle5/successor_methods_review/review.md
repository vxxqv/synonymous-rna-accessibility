# Successor methods and validation-source decision

## Charter

Goal: rank a bounded successor question that can distinguish a measured multi-trait path restriction from additive antagonism, assay transformation and measurement noise. Scope/object/subjects: benign retrospective functional assays / identifiable context or route prediction / paired context transfer, latent-trait nulls, rank tests and count-level uncertainty. Type: exploratory method/source selection, not outcome analysis. Questions: Q1 what new estimand is identifiable; Q2 which priors already occupy it; Q3 which independent material can validate it; Q4 what is the cheapest decisive source gate? Both Cambray arms are now outcome-exposed and cannot validate a successor method.

Sources ranked before collection: original papers/corrections and author methods/deposits (rank 1); institutional metadata/documentation (rank 2); source-code and metadata inspection (rank 3). Search summaries locate primary texts but do not establish methods or novelty. No phenotype files, raw sequencing, fits, installations, external writes or author contact were used.

## Canonical sources

| ID | Primary source / current status | Verification and role |
|---|---|---|
| S1 | Shanas et al., [Promoter identity shapes splicing outcomes and fidelity](https://doi.org/10.1038/s41467-026-76370-1), Nature Communications, 4 August 2026 | Published primary XML, Methods, Discussion; live ArrayExpress metadata, IDF and SDRF read. New context source. |
| S2 | Hidalgo-Carcedo et al., [Allosteric and energetic remodeling of a PDZ domain by protein domain extensions](https://doi.org/10.1038/s41467-026-69673-w), Nature Communications, 19 February 2026 | Published primary XML/Methods/Discussion read; context and assay contract, not phenotype tables. |
| S3 | Martí-Aranda and Lehner, [Seven complete comparative maps of allosteric mutations in a protein family](https://doi.org/10.1038/s41467-026-71005-x), Nature Communications, 21 April 2026 | Published primary XML/Methods/Discussion read; possible separately generated validation material. |
| S4 | Carlson, Andrews and Simons, [Distinguishing direct interactions from global epistasis using rank statistics](https://doi.org/10.1073/pnas.2509444122), PNAS, 2025 | Published primary XML, modeling, limitations and methods read. R&R control and its failure modes. |
| S5 | Faure and Lehner, [MoCHI](https://doi.org/10.1186/s13059-024-03444-y), Genome Biology, 2024 | Primary full text from Amass read; model inputs, multiple phenotypes, preprocessing and fit assumptions checked. |
| S6 | Freudenberg et al., [Lilace](https://doi.org/10.1186/s13059-026-03934-1), Genome Biology, 27 January 2026 | Primary full text from Amass read. Published successor to the 2025 preprint, not a preprint-only method. |
| S7 | Tareen et al., [MAVE-NN](https://doi.org/10.1186/s13059-022-02661-7), Genome Biology, 2022 | Primary publisher model/constraint sections read. Earlier nonlinear-measurement/noise precedent; published version permits scalar continuous responses, not vector responses. |
| S8 | Faure et al., [Mapping the energetic and allosteric landscapes of protein binding domains](https://doi.org/10.1038/s41586-022-04586-4), Nature, 2022 | Primary publisher and institutional author version inspected for assay/data contract. Already establishes abundance/binding biophysical decomposition. |
| S9 | [Mordstein et al. 2020](https://doi.org/10.1016/j.cels.2020.03.001), plus Mikl et al. [2019](https://doi.org/10.1038/s41467-019-12642-3) | Existing authenticated Mordstein audit reused; Mikl2019 exact identifier authenticated through S1 reference 10, not newly reanalyzed. Closest earlier coding/splicing context and library precedents. |

No correction/retraction was encountered in the inspected primary records; the retained XML files have no `related-article` correction link, and the successful Amass records report no retraction. A bounded correction-title search was performed. This is not proof that no correction exists. Amass reached its account limit during the earlier-precedent lookups; two rejected receipts are retained and no further Amass calls were made after observing the rejection. Ordinary public institutional primary XML and metadata access remained available. No restricted access was bypassed.

## Findings ledger

| Finding | Evidence / locator | Consequence |
|---|---|---|
| F1. The existing point path property is ordinal within each trait. | Frozen protocol endpoint minimum/weak inequality and combined interpretation review. | A common strictly increasing transformation cannot remove the same fixed-field obstruction. Fitting an inverse assay curve is not by itself a discriminating biological test. Preparation averaging, filtering, ties and genotype-dependent corrections remain different issues. |
| F2. Additive antagonism already produces joint-only obstruction. | Existing reviewed square example and protocol amendment. | A barrier is not evidence of epistasis. Any successor must add a separate falsifiable model/context question. |
| F3. A new matched-context source is public. | S1 Methods, ArrayExpress E-MTAB-16553 and author repository README. | Same cassette library is measured across promoters; integrated and transient contexts exist. This does not authenticate a synonymous full-ORF graph. |
| F4. S1 already analyzes promoter dependence, cryptic isoforms, predictive sequence features and integrated contexts. | S1 Discussion; Classification Methods uses stratified five-fold CV. | Rediscovering promoter effects, GC dependence or a random-variant prediction score is insufficient. Family/context-heldout calibration is a narrower possible gap, not a priority guarantee. |
| F5. S2 supplies paired functional readouts across four domain boundaries. | S2 Library design, Thermodynamic modelling with MoCHI, Data availability. | A concrete amino-acid method benchmark is feasible in principle, but is outside synonymy and already has context-dependent energetics as prior knowledge. |
| F6. Separate homolog material exists. | S3 Methods and data statement: three biological replicates, raw pre-normalization scores used for modeling. | Potential new material for prediction assessment; same lab/assay lineage and homology still require explicit dependence controls. |
| F7. Rank methods do not solve arbitrary latent dimension or noise. | S4 Model Misspecification and Discussion; single-mutant rank errors, conditional independence, multidimensional GE limitations. | R&R rejection can reflect an omitted latent trait, not direct interactions. Poisson resampling cannot automatically absorb shared batch/normalization error. |
| F8. Multi-phenotype nonlinear models already exist. | S5 model design and Methods; S7 latent measurement process and dataset constraints. | Using a neural network or multiple additive traits is established. The new contribution would need genuinely held-out path/context calibration, not the modeling architecture. |
| F9. FACS-specific uncertainty needs its count data and control contract. | S6 Data processing and Model details; author README. | Lilace can model bin-level overdispersion and replicate baselines, but cannot recover unobserved cross-assay covariance or assert synonymous controls are neutral in a study of synonymous effects. |
| F10. Current measurements cannot validate a redesigned model. | Combined review; both Cambray arms now exposed. | Keep current failed milestones intact; any reuse is explicitly exploratory method development. |

### New source contract, without outcomes

S1 is the most promising RNA-context metadata lead. E-MTAB-16553 API, IDF and SDRF are retrievable. SDRF contains 100 file rows representing **50 distinct source samples**: 40 targeted-reporter samples and ten bulk RNA-seq samples. Paired-end mates are not extra preparations. Promoter conditions have labeled replicate 1/2 samples; Methods authenticates duplicate library cloning. Exact covariance and independence of each integrated preparation are not thereby fully established. The live API reports release 3 July 2026; static IDF says 30 December 2026 despite public access. Preserve this metadata inconsistency.

The reporter cassette is inserted at the GFP coding terminus **before its stop**, so neither untranslated status nor synonymous full-protein identity may be presumed. The reported PSI method counts exon-included reads over all passing reads, with isoform/coverage filters; it is not automatically canonical included/(included+skipped), nor a splicing rate. A canonical-isoform/fidelity endpoint needs a separate authenticated definition. The notebook is 6.56 MB with possible embedded outputs and was **not opened**. No processed CSV, supplement outcome row or raw FASTQ was read. Thus exact construct intersection, synonymy, mutation-edge coverage and family-blocked prediction overlap remain gates, not completed checks.

S2 deposits sequencing as GSE299757 and processed material at [Zenodo 16411462](https://doi.org/10.5281/zenodo.16411462), with [pdzext](https://github.com/lehner-lab/pdzext) code. S3 deposits sequencing as PRJEB90582 and replicate-level fitness/free-energy material at [Zenodo 15978782](https://zenodo.org/records/15978782), with [PDZ_homologs](https://github.com/lehner-lab/PDZ_homologs) code. Their files are named by the primary data statements; no outcome archive was opened. These are measured protein-domain variants in yeast, not same-protein synonymous RNA measurements. Their inferred free energies must not be treated as independent measurements validating the very thermodynamic model used to infer them.

## Scope filter

Discarded unrelated clinical, pathogenic/viral, drug-optimization and population-association discovery hits. No such outcome analysis is proposed. Retained only primary retrospective assay/method evidence. PDZ and Kir2.1 sources below are **method controls**, explicitly outside synonymous-RNA biological validation. No new designs or laboratory procedures are proposed.

## Object filter

Revised the leading default from “run a paired-CDS analysis now” to “resolve a single metadata contract first.” Mordstein's exact-CDS pairing/variant-index gate remains closed; this review found no new evidence resolving it and did not retry the failed supplement. S1 provides a different, accessible context source but cannot silently replace CDS synonymy with a minigene context claim. No method below makes independent assay replication appear from two same-library holdout arms.

## Conclusions: four ranked routes

### 1. Paired context portability, with exact sequence matching first

**Rank / decision:** best next source-only gate; no outcome run yet. The identifiable estimand is held-out context transfer of an explicitly defined assay admissible set or ordered effect, after matching the **same cassette/CDS sequence** and controlling simple composition/splice features. It is not discovery that introns or promoters affect expression. S9 already occupies the GC-by-splicing result; S1 occupies promoter-specific splicing and feature prediction (F3–F4).

**Minimum contract / feasible source:** exact construct index, context-specific transcript start/ORF/isoform annotation, barcode-to-sequence mapping, promoter/preparation IDs, isoform counts/filters, and predeclared family groups. S1's public E-MTAB-16553 and source repositories are a feasible contract-recovery source. Use a strictly synonymous subset only if translation and complete edge/square metadata authenticate it. Otherwise the only permissible endpoint is matched-cassette context portability, requiring explicit scope change. Mordstein remains a reserve rather than an executable substitute.

**Cheapest discriminating test:** metadata-only intersection of exact sequences, reference-exon families, contexts and replicate labels, plus source-code-only review of the original prediction split. Before outcomes, choose one context contrast and retain an independent material/context assessment if overlap permits. Require family-blocked assessment; compare context-invariant, simple GC/length/splice-strength and a training-only monotone context-calibration baseline. A claimed sequence-context interaction must exceed a mere global mapping between contexts.

**Success / failure / stop:** success is reproducible predictive improvement or calibrated nonportability on an untouched declared context/family set, with all missing cases retained. Stop if matching, preparation identity, usable family count, independent assessment or a distinct estimand beyond S1/S9 cannot be authenticated. No externally justified practical margin is supplied by this review. **Cost:** low for metadata, moderate for a bounded count-level analysis; raw remapping is not authorized or required for the first gate. Novelty remains conditional on comparison with S1's full predictor analysis.

### 2. Multi-trait path prediction beyond a specified latent additive model

**Rank / decision:** strongest concrete independent **method-validation** candidate if protein-domain scope is acceptable. S5 already models shared folding and binding traits with nonlinear assay maps; S8 already separates abundance and binding; S2 already establishes extension-dependent energetics/evolvability. The potential gap is calibration of held-out **observed joint path decisions**, not rediscovery of allostery or context dependence (F5–F8).

**Estimand / contract:** incremental out-of-sample accuracy and uncertainty coverage of path acceptability for measured local genotype graphs, relative to a prespecified multi-trait additive/thermodynamic null and a simple local predictor. Require matching nucleotide/amino-acid identities across aPCA/bPCA, measured WT/singles/doubles completing the tested graph, replicate-level scores/counts and errors, explicit boundary context, and training-only preprocessing. Do not validate energies against energies inferred from the same model.

**Feasible source / cheapest test:** first count exact two-assay four-genotype intersections from S2 metadata; reserve genuinely separate S3 material by homolog/ligand and measured family before any outcomes. These are separately generated experiments, not independent laboratories. The full S2 core sequence is not identical across boundary contexts and S3 homologs are not identical proteins; label transfer accordingly. No statement about RNA synonymy follows.

**Success / failure / stop:** the null must allow opposing additive traits and assay nonlinearities. Improvement must occur on held-out measured graph decisions, not only phenotype R-squared, and must survive count/replicate calibration. Stop if completed graphs or overlap are too sparse, batch/source normalizations cannot be reconstructed, or the proposed path statistic is only a reformulation of the published landscape finding. **Cost:** moderate metadata/processed-count work, potentially substantial neural-model fitting later. Model fitting is not authorized here.

### 3. Ordinal falsification of a scalar global-nonlinearity explanation

**Rank / decision:** useful narrow control, not a stand-alone biological accessibility discovery. S4's R&R tests rank constraints from scalar monotone GE; S7 is an earlier explicit measurement/noise framework. The estimand is whether a specified **one-latent-trait** explanation fails, not whether a barrier proves physical interactions (F1, F2, F7).

**Minimum contract / source:** enough repeated mutation backgrounds to estimate rank neighborhoods, input/output counts or a validated error model, precise single-mutant measurements, and material-level replication. Public S4 code and S2/S3 genotype/count material provide concrete method sources; the original S4 datasets are not evidence of synonymy. A handful of selected squares does not meet the large-combinatorial-data premise automatically.

**Cheapest discriminating test:** source-code/method comparison followed, only under a new frozen protocol, by known-answer calibration under additive antagonism, scalar monotone GE, two additive latent traits and shared noise. Rejecting scalar GE must be challenged against the two-latent alternative before any interaction claim. Then assess the fixed statistic on a held-out biological preparation/source with no threshold or rank-neighborhood tuning.

**Success / failure / stop:** success is a correctly scoped scalar-null falsification with reliable error calibration. Stop or narrow the claim if an additional additive latent trait explains the same rejection, single-rank noise dominates, or rank neighborhoods are too sparse. Neither rank invariance nor a deep embedding supplies causal identification. **Cost:** low for analytical controls; moderate for calibrated resampling on adequate count data. The core rank method is published, not novel.

### 4. Count-level uncertainty before joint admissibility

**Rank / decision:** necessary measurement component, not sufficient successor biology. S6 supplies a FACS-specific hierarchical model, whereas DiMSum2020 and Enrich22017 are earlier growth/count-error precedents identified in S4/S6. Lilace models replicate-specific bin baselines and overdispersion; it requires an appropriate control distribution (F9).

**Estimand / contract / source:** calibration of marginal and simultaneous trait/vertex classification errors from bin counts, sort proportions, variant/barcode identities and genuinely traceable preparations. S6's public package and its original Kir2.1 multiparametric source [Coyote-Maestas et al.2022](https://doi.org/10.7554/eLife.76903) are feasible count-model benchmarks. Kir2.1 primarily supplies single-amino-acid effects: it does not automatically supply synonymous two-step squares. Genotype aggregation and barcode linkage must be checked before any graph claim. Do not assume all synonymous variants are null controls when synonymous effects are the scientific target; authentic WT controls are an alternative only if present.

**Cheapest test / criteria:** metadata-only confirmation of bins, sort totals, control identities, independent preparations and an untreated held-out calibration set. Later, compare marginal coverage and simultaneous vertex/path errors under common references rather than independently resampling variants. Calibrated simultaneous marginal bounds can suffice conservatively without paired cross-assay cultures; pairing is not logically indispensable. However, neither uncalibrated marginal intervals nor separate posterior draws magically identify joint biological covariance.

**Stop / cost:** stop a path claim if there are no measured intermediates or no trustworthy controls/uncertainty calibration. Do not bootstrap the existing processed Cambray columns into fictitious biological replication. Cost is low for schema review, moderate-to-high for a carefully validated hierarchical fit. Use this component with routes 1–3 only after their source gates pass.

## What the current methods cannot solve

- A transformer, protein/RNA language model, active-learning acquisition score or larger neural network cannot identify unseen measurements, restore a masked target's outcome or create independent preparations. Predictions are not validation observations.
- A fitted monotone assay correction does not remove a fixed-field ordinal barrier. A nonmonotone map can change it but requires separately justified assay physics, not favorable fitting.
- R&R rejection does not identify direct interactions in a multi-latent system. MoCHI's good fit does not establish the correct latent dimension or causal biological mechanism.
- Marginal Bayesian credible intervals do not automatically calibrate a simultaneous multi-vertex path certificate or shared assay effects. Borrowing across positions/variants is a model assumption to validate.
- A random variant split can leak shared exon families, sequence backgrounds, normalization and targets. A context holdout in one library is not an independent experiment merely because its values were unread.
- More single mutants, convenient thresholds, a prevalence estimate after filtering, or a new name for an existing context effect do not establish multi-step accessibility.

## Next steps / Risks

**Recommended next decision:** authorize at most the route-1 metadata/source-code gate, with no phenotype values. Recover exact construct/annotation and sample contracts, compare the original predictor's split, and decide whether a genuinely different paired-context question survives. If strict synonymous coding identity is indispensable and fails there, stop that branch; route 2 is a separately labeled general-method benchmark, not a substitute biological validation. No route presently clears every independent-validation, calibration, same-CDS and graph-coverage gate.

The review is bounded, not exhaustive: three Amass discovery lanes, targeted primary/correction/author-source searches and directly identified public institutional texts/metadata. Search and exposure records, failed routes, file hashes and source-specific limits are in the companion records. The current negative/information-limited Cambray conclusion remains unchanged.
