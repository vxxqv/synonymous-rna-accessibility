# FAS synonymous accessibility pilot: protocol draft

Status: DRAFT. Outcome access is prohibited until this document, the analysis implementation, synthetic verification and independent source review have been frozen together. This is an adaptively selected retrospective analysis of published data. It is not a preregistered experiment or untouched independent validation.

## Question and scope

Does the measured enrichment ordering within the entire declared native-reference-protein FAS subset require more than arbitrary codon effects followed by a common monotone response? If so, does that departure coincide with a reproducible accessibility restriction inside the declared graph?

The endpoint is relative included-RNA enrichment, a source assay proxy for exon inclusion. It is not calibrated absolute PSI, endogenous transcript function, measured protein abundance, apoptosis, organismal fitness or evolutionary selection. Variant-dependent RNA abundance, stability, RT or PCR recovery can alter this endpoint and are not removed by normalization alone. Native included-isoform protein identity is authenticated from NM_000043.6 and NP_000034.1; actual translation of the experimental minigene and preservation of its isoform mixture are not claimed.

## Fixed domain and source

Use the full 72 sequences enumerated before outcomes from the published 3072-sequence design at Scaling_Law commit b87064c39a553b0759d86c2b428813c58dc17ae2. The fixed genotype list is the native-reference translation class, not the ancestral translation class. Its graph has single-nucleotide edges only. The expected product graph is P3 x P3 x Q3, with 204 edges. Independent enumeration must confirm these facts before outcome access.

Read only the twelve pinned count files documented by the source: nine biological output transfections, three technical DNA-input PCRs. Preserve complete source-file totals and the author's TR1-aligned totals for reproduction of the original median-input normalization. The author aligns rows to TR1 and omits missing values from sums and medians; explicitly record any difference between those totals and all-file totals. Zero aligned totals make those normalized frequencies undefined and are recorded; they do not abort the primary missing-data analysis. Do not normalize only within the 72-sequence subset. Do not run raw-read reprocessing or download new samples. No new source may replace this domain after outcomes.

The source processing README exposed five illustrative genotype/count rows. They are recorded in the exposure ledger, cannot be described as untouched, and will not determine selection or endpoints. Published findings, source Methods and aggregate coverage statements are also exposed.

## Measurement and uncertainty

Let O[v,r] and I[v,t] denote deposited output and input counts. Let w be the human reference. For every output r and input t with positive required denominators, form

R[v,r,t] = (O[v,r] / O[w,r]) / (I[v,t] / I[w,t]).

This gives 27 normalization scenarios sharing preparations. They are not 27 independent biological replicates. The original author's per-output ES, computed from normalized output frequencies divided by median normalized input frequencies, will be reproduced separately and then divided by the reference ES in the same output. The descriptive point summary is the arithmetic mean of these nine relative ES values, provided all nine are defined. It is not the ratio of the two mean absolute ES values or an averaged calibrated PSI. No numerical PSI calibration or author's outcome-selected filtering is borrowed.

For each v, the primary envelope [L[v],U[v]] is the minimum and maximum over these scenarios. Finite known zero outputs are allowed and yield zero enrichment. No pseudocount is used in the primary analysis. Any absent record, unauthenticated duplicate, missing value or zero input denominator makes that vertex unknown with envelope [0,infinity] unless a stricter rule is authenticated and frozen before access. A nonpositive reference count makes the affected scenario undefined and the primary full-scenario envelope unknown, rather than silently dropping that scenario. Malformed/negative/noninteger counts or conflicting duplicate records stop the pilot.

These envelopes describe observed preparation/normalization variation. They are not confidence intervals, do not bound unknown RT/PCR bias, and do not establish sampling coverage for a population mean. The same library, shared inputs and sequencing lane limit independence. No p-value, population-level confidence statement or biological-equivalence claim will be generated from the envelopes. Source-supported calibration would require a separate protocol before any inferential claim.

## Primary null and decision

The null is R[v] = h(X[v] b), where X contains an intercept and categorical effects for each of the five changing codons, and h is any common nondecreasing function. X has eight linearly independent columns on the complete design. The coefficients are not uniquely identifiable under an unrestricted h, including location and scale freedoms. Arbitrary within-codon effects on the two three-state paths are permitted. A dip at the middle state of one codon can disconnect this graph under the null.

Construct an order constraint (X[j]-X[i]) b >= 1 for each exact forced order U[i] < L[j]. With nonempty closed envelopes, feasibility is equivalent to some common nondecreasing response taking all predictions inside the envelopes. The finite strict constraints can be scaled to margin one; this margin is a mathematical normalization, not a measured biological effect size. Infinite upper bounds impose no outgoing forced order and must be handled explicitly.

If feasible, return an inspectable rational vector, check every margin exactly, and verify an assigned common monotone response inside every envelope, including latent ties. Rationalization may be rescaled by its exact positive minimum margin; otherwise the result is unresolved. State that the null has not been excluded, not that it predicts new genotypes. If infeasible, obtain a nonnegative dual certificate with A-transpose times lambda = 0 and sum(lambda)=1. Convert the certificate to rational values and check all equalities exactly before reporting a certified contradiction. Floating optimizer status alone is insufficient. Failure to verify a certificate is an implementation/information limitation, not a scientific rejection.

Primary result: feasible, exactly certified infeasible, or unresolved because the solver/measurement contract failed. A missing-data graph can remain informative, but no count of missing genotypes is relabeled as a biological negative. No outcome-selected subset or fitted sigmoid replaces this null.

This is an application of established finite order representation and linear inequalities. Closest-prior comparison must include the FAS global-scaling analysis, monotone global-epistasis models, ranking approaches and additive conjoint measurement. No new foundational theorem or folding-specific mechanism is claimed.

## Accessibility and effect size

For every unordered endpoint pair in the full declared graph, compute the maximum bottleneck score over paths, using lower and upper envelope fields separately. A conservative obstruction requires endpoints certified above a threshold by their lower scores, while even the upper-score graph cannot connect them at that threshold. With acceptance defined by score >= threshold, the excluded-connection threshold set is (upper bottleneck, minimum endpoint lower score]. Equality of its two endpoints gives an EMPTY set; an infinite upper bottleneck also gives no certificate. Lower bottleneck >= threshold guarantees a connection within the assumed envelopes. All other cases remain unresolved. Report all pairs, including empty threshold intervals. No favorable threshold or endpoint is selected as the sole report.

Also compute the same descriptive quantities for every defined single normalization scenario and the reproduced author point estimate. Preserve all results. An obstruction applies only inside the 72-state design. It does not exclude synonymous routes through unmeasured alleles, positions or higher-dimensional contexts.

Report witness contrasts and ratios on their original scale. No arbitrary fold cutoff becomes a claim of biological importance. A positive order contradiction without an accessibility restriction supports only the narrower ordering result; a barrier without a contradiction may be explained by codon gates. Even the co-occurrence of both does not show that the particular barrier requires intercodon interactions: one codon can cause the barrier while another part of the graph violates the null. No barrier-specific mechanistic attribution is allowed without an additional targeted necessity certificate. Both outcomes, including neither, remain reportable.

## Coverage and model checks

Retain all 72 genotypes. Report counts of finite envelopes, unbounded envelopes, observed records and unknown states separately. Expected graph membership is fixed by sequence design and must not be reduced to observed vertices. Full measurement means all 72 vertices are identified; otherwise only explicitly bounded/partial results may be reported. No minimum sample number is changed after observing coverage.

On this particular 72-state design, nucleotide-additive and arbitrary codon-unary latent fields have exactly the same column space. Each three-state codon factor occupies three corners of two binary nucleotide coordinates, so its intercept plus two nucleotide effects spans all three codon values. Verify the rank-eight column-space equality exactly before access. The two encodings are an implementation equivalence check, not competing biological baselines. A three-state dip is possible even under nucleotide additivity on this restricted domain; transferring full binary-cube connectivity to it would be wrong. The full 72-row raw additive design has 64 independent linear contrasts. These are not statistical residual degrees of freedom for the unknown-monotone-response model, nor independent biological units; measured coverage is not yet known. A nonlinear h can generate nonzero raw additive residuals even under the null. Include the saturated genotype model only as an adequacy reference, never as predictive generalization. Do not fit extra interactions or machine-learning models to rescue a failed result. Subsequent predictive work would require a separate frozen split and independent biological material.

## Required pre-outcome synthetic checks

Verify: a binary additive field and arbitrary nondecreasing response are accepted; a P3 middle-state dip can create a barrier yet is accepted by codon additivity; a two-codon reciprocal-sign fixture is rejected with an exact dual certificate; constant/tied fields are accepted; overlapping intervals can render a point contradiction unresolved; an unknown bridging vertex cannot produce a false full-graph exclusion; row permutation and a common strictly increasing transformation preserve order feasibility; duplicate features with disjoint intervals are rejected; empty forced-order systems produce a valid constant response. A merely nondecreasing transformation may collapse orders and remove a contradiction. Compare bottleneck outputs with an independent exhaustive/reference calculation on small graphs. Check the shared-input normalization separately.

## Validation and stopping

The nine outputs are repeated transfections of one shared library. This pilot uses all nine for conservative descriptive envelopes; none is labeled an independent validation cohort. FAS2016 can later test limited same-exon reproducibility, but its overlap, design influence and normalization must be frozen before its outcomes. It does not supply another independently sampled exon.

Stop after one full-domain analysis and independent computational review. Do not retune envelopes, thresholds, codon factors, genotype membership or source selection to obtain a positive result. A feasible null, no barrier, insufficient measurements, and uncalibrated biological uncertainty have distinct interpretations and must remain distinct in the decision record.
