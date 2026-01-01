# CNOT3 association, baseline stability and depletion response

The focused investigation is complete as an explicitly exploratory prediction comparison. Adding measured CNOT3 association produced a small improvement in the recovered population, but failed the frozen practical performance criterion. The result concerns deposited lifetime estimates within one study. It does not reject the established CNOT3 mechanism or identify a pathway-specific decay rate.

## Question and design

Does measured association add prediction of depletion lifetime beyond baseline lifetime, codon or amino-acid composition, transcript features and available measurement covariates? The primary endpoint was natural-log depletion lifetime. All preprocessing, projection of association onto the codon comparator, ridge tuning and prediction used training genes only. Near-protein groups, rather than individual genes, defined the fixed split.

The published 5,481-gene population was not reproduced exactly. The deposited table supplies 5,480 positive paired half-lives, and the additional historical identifier remains unresolved. Forty-nine stored log2 changes fail the original arithmetic tolerance. A post hoc precision diagnostic reproduces all stored changes from ten-significant-digit half-lives, but the undocumented generating operation is not proved. The original reproduction gate remains failed. No gene or value was corrected to repair it.

One operational-cohort successor was recorded after that source diagnosis. It retained the endpoint, exclusions, model definitions, partition salts, tuning, sensitivities and 5% performance criterion. Source reproduction had already exposed the cohort outcomes. The result is therefore model-held-out prediction after source exposure, not prospective validation or independent experimental replication. The later reanalysis and the original CNOT3 study use the same underlying biological material.

Reference eligibility retained 5,416 of 5,480 genes (98.832%). The fixed split contains 4,304 training genes in 4,280 protein groups and 1,112 held-out genes in 1,109 groups. An independent cohort audit regenerated every candidate relationship, exclusion, predictor and partition. The 72 near-protein candidate distances were checked exactly. This protects the declared sequence boundary, not dependence across experiments, pathways or regulatory states.

## Complete prediction results

| Model | Held-out MAE | Held-out RMSE | Spearman correlation |
|---|---:|---:|---:|
| No change | 0.473567 | 0.633015 | 0.806023 |
| Linear baseline | 0.362495 | 0.484391 | 0.806023 |
| Nonlinear baseline and core | 0.341733 | 0.456460 | 0.829597 |
| Amino-acid and core | 0.339458 | 0.455163 | 0.830842 |
| Codon and core | 0.339700 | 0.455042 | 0.831840 |
| Association, primary | 0.339316 | 0.455127 | 0.832215 |
| Codon without enrichment SE | 0.339756 | 0.455048 | 0.831871 |
| Association without enrichment SE | 0.339418 | 0.455141 | 0.832195 |
| Fixed nonlinear association | 0.338615 | 0.454404 | 0.832647 |

MAE and RMSE are in natural-log lifetime units. The primary gain is 0.113113% over codon composition, an absolute MAE reduction of 0.000384244. It is only 0.041944% over the strongest primary baseline, the amino-acid model. The criterion of at least 5% improvement over every primary baseline fails. Primary RMSE is slightly worse than the codon comparator.

The two fixed sensitivity gains over their matched codon comparator are 0.099518% without enrichment SE and 0.319333% with nonlinear association. They do not replace the primary result. All nine prediction vectors and all 10,008 held-out model/gene error rows are retained.

Deleting each of the 1,109 held-out protein groups without refitting gives codon-relative primary gains from 0.078287% to 0.139775%. The gain remains positive against codon in every deletion, but not against every baseline. These are finite-cohort influence diagnostics, not confidence intervals, biological-replicate uncertainty or proof of transportability. No calibrated half-life uncertainty was available.

## What baseline lifetime and rate changes identify

The unit-invariant reciprocal contrast is R = ln(2) t_ref (1/t_control - 1/t_depletion), with t_ref the training median control lifetime, 5.17814182006782 in the deposited units. Across all 5,416 eligible genes, its Spearman correlation with log2 lifetime change is 0.826797. The respective correlations with baseline lifetime are -0.581621 and -0.141675. These algebraically coupled endpoints therefore rank genes differently. Common rescaling of all lifetimes leaves R unchanged; this is not a calibrated molecular rate.

All 912 negative reciprocal contrasts remain in the data. They are not 912 demonstrated compensatory effects or negative pathway contributions, because point-estimate error and processing conventions remain unresolved.

Under single-exponential kinetics, k0 = b + p and k1 = b + delta + r p separate a baseline pathway component p, background b, residual activity r and background change delta. The net rate reduction D = k0 - k1 equals (1-r)p - delta. Without justified bounds on those nuisance quantities, it does not identify p.

For independent external ranges r in [rL,rU] and delta in [dL,dU], the sharp feasible conditions are 0 <= p <= k0, (1-rL)p >= D+dL, and (1-rU)p <= D+dU. A nonnegative perturbed background additionally requires rL p <= k1 and p <= k0+dU. If zero residual activity remains admissible and delta is unrestricted, every p in [0,k0] is compatible with the same positive observed rates. This last statement must not be asserted under arbitrary positive lower bounds on residual activity. These are elementary conditional identified sets, not new general kinetic theory or confidence intervals. No source-specific nuisance bounds were recovered.

## Verification and research decision

The operational amendment and independent method review were committed as `9d0691960f81e777eea6de889bdc8229db7b4c6e`. All nine predictions were committed as `b709472c473e1f7ff53a20d8c7fe28ec97b43bf1` before scoring. Prediction SHA-256 is `41661565bf1eb49398fa5cc9141811db4b0f81b10ff77f55a075984897718b93`. Evaluation SHA-256 is `95fb1c15cfe6148c41a543240e1c84db0d54b910cbe974eeec26914202ceb2ac`.

Independent numerical reconstruction passed 59,268 checks comparing 240,144 numerical values, with maximum absolute difference 3.3973e-13. It checks the stated algorithms and deposited calculations; it does not authenticate biological validity or historical access order. The source discrepancy, omitted outcome-error model and unresolved preparation identities remain explicit.

The focused prediction question has an interpretable answer: measured association supplies little incremental accuracy in these fixed models and this operational cohort. Preserve the practical failure, retire the examined holdout from tuning and stop expansion on this dataset. Distinguishing a small true contribution from attenuation due to noisy association, endpoint error or baseline adjustment requires new measurement information or genuinely independent material. More flexible fitting on the examined cohort would not supply that evidence.

## Claim to evidence

| Status | Claim | Supporting evidence | Remaining limit |
|---|---|---|---|
| Established calculation | Association adds only a small practical prediction gain in the declared comparison | All nine committed predictions; 1,112 held-out genes; independent reconstruction | Exploratory after source exposure, one study, error calibration absent |
| Established conditional mathematics | Two decay rates do not uniquely allocate a pathway contribution | Constructive counterexamples and sharp nuisance-bound inequalities | Requires exponential and additive-rate assumptions; no empirical nuisance bounds |
| Exploratory description | Lifetime change and reciprocal-rate contrast rank genes differently | All 5,416 genes, rho 0.826797, all 912 negatives retained | Shared baseline measurements induce coupling; no molecular-rate calibration |
| Unresolved source limitation | Published cohort not reproduced exactly | Independent readers recover 5,480 instead of 5,481 | Missing historical identifier and original processing unavailable |
| Incomplete inference | No causal pathway allocation, folding mechanism or genotype accessibility conclusion | Eligibility and identifiability audits | Needs appropriate measured endpoints, error model and independent biological design |

The complete tables accompany this report as CNOT3_All_Model_Results.csv, CNOT3_Verified_Heldout_Predictions.csv, CNOT3_Component_Influence.csv and CNOT3_Kinetic_Comparisons.csv. Original source and arithmetic evidence is retained in CNOT3_Source_and_Kinetic_Audit.md and CNOT3_Arithmetic_Diagnostics.csv.
