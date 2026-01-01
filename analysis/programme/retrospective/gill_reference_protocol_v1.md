# Gill reference-protein pilot protocol

Status: frozen before any RNA outcome is decoded. The six-library source mapping is verified; biological independence and pairing remain unresolved. Independent numerical and selective-exposure review passed at the recorded source hashes. All 27 pilot tests and the complete 77-test project suite passed. Every prediction must additionally be committed and verified before the held-out RNA library is decoded. A failed result is retained without retuning.

## Question and inferential scope

Does a fixed, partially pooled pairwise-interaction predictor reproduce the RNA measurement of a held-out library better than simple sequence features, an additive model, a two-stage monotone additive fit, graph smoothing and the observed training mean on the same 32 genotypes?

The estimand is a reference-normalized contrast of pseudocount-adjusted RNA and DNA read counts from an existing deposited library. It measures relative RNA abundance in this assay. It is not a protein-output, fitness, structural-occupancy or decay-rate endpoint. This is a test of held-out RNA-library reproduction conditional on an already inspected DNA denominator. It is not validation on new genotypes, proteins or experimental settings. The metadata audit maps all six libraries to separate archived samples and experiments, with distinct raw-file checksums. All are labelled pooled, but neither pooling semantics nor independent cultures, extraction identities or RNA/DNA pairing are established. Paired-end sequencing does not imply biological pairing. Unknown dependence will not be replaced by a convenient independence assumption.

## Fixed design and exposure

The eligible reference protein is the deposited 313-residue reporter. Its complete synonymous fibre has 32 existing genotypes, 112 one-nucleotide edges and graph K2 Cartesian K4 Cartesian K4. Authenticated construct metadata, rather than RNA measurements, determine the literal ordered keys and reference key in the accompanying JSON. Exactly those keys are retained. All three DNA libraries were inspected during eligibility; no individual RNA outcomes have been opened before this protocol. Published aggregate findings and source methods are already known.

The fixed source is Mendeley Data version 1, DOI 10.17632/ds7zymk875.1. The processed workbook SHA256 is 0f80ba9868c54548db938c65d09d8caa0419c5721bfa379949f6a18ad2568e89. Construct, parser, amplicon and eligibility hashes are retained separately. Source counts, not historical shrunken fold-change or error columns, are used. The source's DNA coverage selection is part of the conditioning set. Raw-read processing is not independently reproduced in this pilot.

## Training, normalization and holdout

Use RNA libraries CDS1 and CDS2 and DNA libraries CDS8 and CDS9 for training. Use CDS3 RNA with CDS10 DNA for validation. RNA and DNA are not assumed paired biological preparations. CDS10 DNA is already exposed. At pseudocount c, the unanchored training vector is the mean of log2(RNA+c) across the two RNA libraries minus the mean of log2(DNA+c) across the two DNA libraries. The validation vector is log2(CDS3+c) minus log2(CDS10+c). Subtract the value at the fixed reference genotype from each complete vector. Primary c=0.5; c=0.25 and c=1 are frozen sensitivity analyses. Save all three sets of predictions before examining CDS3.

The within-field reference subtraction removes a common library-depth offset after the count pseudocount. It does not remove measurement error in the reference genotype. That error produces shared dependence among the 31 nonreference contrasts. Do not normalize counts before adding the pseudocount, estimate normalization factors from the validation RNA vector, or fit any transformation on that vector.

## Frozen predictors

Fit with equal weights on all 32 training genotypes, then subtract each predictor's own fitted reference value. The full factorial decomposition has dimensions 1, 7, 15 and 9 for constant, main, pairwise and three-factor effects. Nucleotide and codon additive spaces coincide here and are one comparator.

1. Constant, giving the zero contrast predictor.
2. Ordinary least squares on intercept and total GC count.
3. Uniform orthogonal projection onto intercept and categorical main effects, rank eight before anchoring.
4. The same additive predictions followed by weighted least-squares nondecreasing isotonic regression of training response on those predictions. Group exactly equal input predictions, weight by group size, use deterministic PAVA and merge equal fitted blocks. Interpolate at input knots and clip outside their range. This two-stage fit is not an optimization over every monotone-additive representation. It may have up to 32 fitted blocks and is not guaranteed unsaturated.
5. The primary pairwise predictor: keep constant and main effects, retain one half of each orthogonal pairwise coefficient and omit the nine three-factor directions. Under uniform empirical orthonormalization this is ridge penalty one for pairwise terms, with smoother trace 15.5 before anchoring and 14.5 afterward. These are algebraic smoother traces, not biological error degrees of freedom.
6. Graph heat smoothing exp(-tL) of the training vector, with unnormalized Laplacian L and t=log(2)/4, followed by anchoring. Its trace before anchoring is 10.6694173824.
7. The per-genotype training vector itself, a strong saturated measurement-reproduction comparator.

There is no tuning, model selection or omission of a comparator after validation. Every predictor, its settings and all 32 values are retained for every pseudocount.

## Endpoints, hurdle and uncertainty

The primary endpoint is RMSE on all 32 anchored validation contrasts at c=0.5. The reference contributes zero error mechanically. The accompanying 31-genotype RMSE differs by the fixed factor sqrt(32/31) and cannot change model rankings. Report every vertex error, signed mean error and centered RMSE, where the error vector's own mean is removed. Centered RMSE diagnoses common-reference sensitivity and cannot replace the primary endpoint after inspection.

The predeclared descriptive success hurdle is pairwise-predictor RMSE at least 10% below the minimum RMSE across all six non-pairwise comparators, using unrounded values. The best comparator loss must be strictly positive. A zero-error comparator cannot be improved by 10%. Both sensitivity pseudocounts must retain the direction of improvement before the result is called insensitive to the pseudocount. Failure does not license another endpoint, parameter or selected subset.

The 10% hurdle is a practical decision rule, not a significance test. No calibrated sampling model or simultaneous confidence region is established. Do not bootstrap genotypes as independent biological units, use the factorial residual dimensions to estimate biological sampling error, or attach biological p-values to interaction energy. Descriptive prediction errors are exact summaries of this finite held-out measurement vector. Persistent sequence-dependent assay bias can reproduce across libraries.

## Secondary geometry and sensitivity

After primary validation, retain the orthogonal interaction energies, every critical-threshold component count and the maximum endpoint-relative barrier across all 496 distinct unordered endpoint pairs. Repeat the fixed summaries for the nine RNA-by-DNA library combinations. These nine fields share libraries, reuse training measurements and are dependent sensitivity scenarios. Their finite extrema are neither confidence bounds nor guarantees over their convex hull or every future experiment.

Any exact additive field or common monotone transformation of it has connected nonempty superlevel sets on this product of cliques. A positive observed barrier excludes that representation of the observed vector but can be generated by measurement noise. Nonzero interaction energy alone can also arise from a common nonlinear monotone response. Neither observation establishes a folding mechanism or biological interaction.

## Preflight, stop rules and audit trail

Before observed fitting, independent synthetic checks must recover known additive and nonlinear-monotone negative controls and a known interacting positive control, verify the factorial dimensions, and compare graph calculations with a separate explicit algorithm. Review the selective workbook reader to ensure training never decodes CDS3 outcomes or unrelated protein outcomes.

Stop on changed source hashes, changed metadata mapping, headers differing from the frozen schema, anything other than the fixed 32 unique keys, a missing or duplicate reference, missing or nonfinite counts, negative or noninteger counts, or evidence that nominally separate libraries are duplicates or derived copies. Exactly duplicated selected RNA vectors trigger a conservative identity hold. Zero RNA counts alone are valid and are retained. No outcome-based exclusions are permitted.

The primary agent records the frozen protocol and implementation revision before training, saves and hashes all predictions, and commits them before opening validation RNA. The validation record states the exact files and columns first exposed. If the metadata do not establish independent biological preparations, the conclusion remains descriptive library reproduction regardless of performance. Stronger biological claims require additional evidence and are not authorized by a favorable score.
