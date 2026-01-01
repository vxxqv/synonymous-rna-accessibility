# COMT A1 interpretation limits, fixed before outcomes

Status: PRE-OUTCOME interpretation memo. The exact protocol v2 and integration-addendum bytes are pinned in receipt.json. This is not authentication of a final freeze, outcome authorization or an independent biological replication. No COMT count body, source simulated-training value, phenotype, prediction or observed result was opened for this memo.

Let d, r, h and f denote the declared central natural-log representations for gDNA, total RNA, heavy polysome and the low-fluorescence bin. The frozen contrasts are R=r-d, L=h-r and T=d-f. Thus, including the fitted intercept c,

    c+aR+bL = c+(a-b)r+b h-a d.

The constrained joint model has a>=0 and unrestricted b. Its two covariates are derived from three raw assay summaries, not independent biological measurements. In particular, a=b eliminates measured RNA r completely: the fitted variable part becomes a(h-d). A joint-model improvement can reflect a more useful combination of heavy representation and gDNA, cancellation of RNA measurement error, source preparation/normalization structure or true biological differences. These alternatives are not separated by held-out codon sites. No sign or magnitude of b alone identifies a ribosome-mediated causal effect, translation flux or an independent biological load mechanism.

## Existing exact falsifier

The pre-outcome, independently authored control in ../comt_independent_review/shared_reference_checks.py and its saved JSON fixes true RNA=heavy=u, so true heavy/RNA log load is identically zero. Measured RNA is r=u+e, h=u, d=0 and f=-u, with u in {-1,1} and e in {-2,2}. Therefore R=u+e, L=-e and T=u. The joint fit a=b=1 has MAE 0. The four ordered RNA inputs are (-3,-1,1,3), with targets (-1,1,-1,1); nondecreasing least-squares PAVA gives (-1,0,0,1), MAE 1/2. The constant baseline MAE is 1. Replicating these measurement patterns at distinct training and assessment sites would preserve this predictive advantage without any variation in true load. This is a logical counterexample to mechanistic identification, not a proposed explanation selected from observed results.

Even when raw measurement errors are centered and mutually independent, the contrasts share errors: Cov(R,T)=-Var(d) and Cov(R,L)=-Var(r). Biological effects and additional dependence can modify these covariances; these identities describe the isolated independent-error component, not the observed total covariance. A held-out site split cannot erase repeated assay-reference algebra or turn preparation summaries into independent validation experiments.

## Prospective outcome interpretation

| Frozen result pattern | Defensible interpretation | Excluded interpretation |
|---|---|---|
| Coverage >=100 variants and >=40 sites, and heavy-model site-balanced MAE reduction >=10% against each of all four baselines under EACH code/GEO map | The named joint predictor passes the operational expansion rule for this source-ascertained construct and raw-proxy reconstruction on held-out codon sites. Report all errors, preparation summaries, ROI and fixed-SNV diagnostics. | Independent replication, a causal load mechanism, physical protein abundance, assay equivalence, calibrated uncertainty, or superiority over every RNA-only predictor. |
| One map passes and the other fails a magnitude, coverage or stability requirement | Primary operational rule fails; the descriptive comparison is source-map sensitive. Preserve both maps and results. | Choosing the favorable map, treating map agreement as authenticated sample pairing, or treating the two maps as independent studies. |
| Neither map passes, a zero baseline MAE occurs, or common coverage is insufficient | The predeclared pipeline did not meet its expansion criterion. Distinguish numerical/source failure, insufficient coverage and a completed predictive comparison. | No possible RNA predictor can work, load is biologically irrelevant, assays are interchangeable, or missing measurements represent biological zero. |
| Mid-polysome or fixed 113-SNV sensitivity is favorable after the primary fails | A declared descriptive sensitivity, reported with the same restrictions. SNV errors use the same predictions without refitting. | Rescue of the primary, a new selected endpoint/cohort, or confirmation of a mechanism. |

The target remains negative low-bin enrichment, not mean fluorescence or protein concentration. Both maps reuse the same underlying source measurements with different preparation-group labels. Their agreement is sensitivity robustness to these two named maps, not proof either map is correct or exhaustive. The fixed source preprocessing, RF reconstruction and published-score ascertainment remain conditioned inputs. Site cross-validation assesses calibration across source sites, not transfer to new preparations, constructs, native isoforms or populations.

The four baseline requirement does not make the comparison exhaustive. PAVA minimizes training squared error but is assessed by MAE; its interpolation/clamping differs from linear extrapolation. The RNA-linear model is nested in the joint model for training squared error, but held-out MAE improvement is not guaranteed. The 10% and coverage thresholds are operational decisions, not biological equivalence margins. Preparation envelopes are descriptive sensitivity sets, not confidence intervals or protection against persistent source bias. No outcome-dependent model, threshold or alternative is introduced here.

This memo fixes interpretation of either outcome before any authorized COMT numerical reconstruction. It deliberately makes no empirical prediction about which alternative will occur.
