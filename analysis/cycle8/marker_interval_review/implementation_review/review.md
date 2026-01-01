# cPRIME B1 pre-outcome implementation review

Verdict: PASS for the current pinned, dataset-specific implementation and descriptive protocol, conditional on the root completing an immutable freeze of these code/protocol/eligibility bytes and all declared source hashes before authorizing any real outcome access. No real outcome access or production run was performed in this review. No source files were edited by the reviewer.

## Resolved findings before freeze

| Finding | Prior review point | Current resolution |
| --- | --- | --- |
| Finite nonpositive source scores skipped reproduction comparison despite positive paired RPM | runner b7fa03dc8da8440fcaa55139d566573c49dc6b32bf62bbffee000010400f529a | runner first revised to 540f5c4ff7fb6f6ca2ec2145e0b6b0137ece1a7dbe6fb915d1029851347aa9ae, then retained in final bfba04ed version; independent zero/negative score fixtures now return SOURCE_REPRODUCTION_FAILED |
| Missing measurement versus absent metadata row was ambiguous | protocol b709f46bdab736bdb4d50ca9948ea4b986ef3378b751e4bbff6deea8945ac145 | current protocol distinguishes unknown measurement interval from integrity-fatal missing/extra sealed metadata row |
| Scientific module imported before authorization/hash verification | intermediate runner 540f5c4ff7fb6f6ca2ec2145e0b6b0137ece1a7dbe6fb915d1029851347aa9ae | deferred until after main gates; modified-module import canary remains unexecuted in isolated test |
| All shared-string text and header fields eagerly materialized | same intermediate runner | only indexed shared strings and expected headers are semantically decoded; protected shared-string tracer and unselected formula fixtures pass |
| Outcome conversion preceded validation of eligible row context/identity | same intermediate runner | identity is now rejected before conversion of poisoned outcome values |

Original observed hashes and findings are retained here rather than treating earlier revisions as accepted. The root also reported a prior None-annotation preparation failure and a Windows synthetic-file-handle failure corrected before values; this review does not independently reconstruct those historical executions. They are not empirical analysis failures or outcome-driven exclusions.

## Current mathematical/protocol checks

For positive response ratios, the minimum multiplicative spread is max(1,max L/min U), and the supremum uses distinct markers max(U_i/L_j). A single marker has spread 1. Unknown (0,infinity) intervals add no forced order and, for a multi-marker set, make the upper spread unbounded. Finite positive bounds are attained; infinity is an unbounded supremum, never an observed response. The implementation handles partial zero lower bounds consistently with strict positivity. All order/claim decisions use rational arithmetic; no logarithm or floating tolerance changes these decisions.

The unrestricted scalar-order null is exactly the protocol's stated common latent score with context/exon-specific nondecreasing maps. Every union edge is supported by a strict rational upper/lower separation in an identified cell. A directed cycle contradicts every such latent order; an acyclic graph permits a monotone completion. The small positive starting value in the returned completion is below all positive finite bounds, so it handles all-unknown and zero-lower-bound cases without inventing a zero response. A cycle is not an interaction, mechanism, additive-model fit or causal context effect.

Independent completion enumeration passed 2,379 patterns from 13 positive-domain intervals, including fractional endpoints and unknown/partially unbounded intervals. A separate value-assignment/permutation oracle passed 2,048 two-context cases, checking every returned edge origin, all feasible response witnesses and all cycle edges. This is supplementary finite evidence; the general argument is the elementary interval/order proof, not the finite enumeration.

Metadata checks confirm exactly 30 groups, 66 states, 24 two-marker and six three-marker groups; group memberships partition the eligible states. The eligible exons are 17,25,27,39,42, not six exons. All currently eligible state annotations are exactly synonymous_variant. The preparation script is dataset-specific and its substring annotation filter should not be generalized to arbitrary future annotation vocabularies without another review.

## Reader and pipeline challenge

The independent invented full-size fixture uses the eligible metadata identities with synthetic values only: 396 rows, 132 state/cell records, 60 comparisons and five exon tests at each of four widening factors. All exact bounds were independently reconstructed by completion enumeration; all five-exon statuses were checked through a separate transitive-closure cycle calculation. Row permutation preserved scientific results. Within-row paired RNA/DNA division was checked using ratios 4,1,1/4, rather than a ratio of separately pooled or permuted quantities. No additional pseudocount is applied to deposited RPMs.

The reader preserves all required rows, marks blank/zero/nonfinite required measurements as an unknown interval and retains all states. Missing or extra sealed metadata rows fail integrity. A finite zero or negative score is compared to a valid paired ratio before usability filtering. The exact equality boundary of the 1e-8*max(1,abs(q)) consistency test passes and an immediately larger rational disagreement fails at both q=1 and q=100. This tolerance is a file-reproduction check only.

Protected unselected body rows, unused outcome/formula columns, an unused header formula and an unselected shared-string text tracer were not semantically decoded. A wrong-context eligible ID was rejected before a poisoned numeric outcome was converted. XML byte parsing is not claimed to avoid reading the containing file bytes; the tested restriction is semantic selection and numeric conversion.

Only an isolated copied CLI layout with an invented 396-row workbook was executed. Matching toy authorization/freeze/source hashes produced the expected complete result. Denied authorization, wrong freeze hash, changed scientific-module bytes, changed synthetic source bytes, existing output and optimized Python (-O/-OO) were rejected. The modified-module canary did not execute, demonstrating the corrected import order. The original source workbook was never opened or hashed by this review script.

## Scientific claim boundary and remaining freeze conditions

The source-derived endpoint is a local motif RNA/DNA proxy. Empirical preparation envelopes are a specified Cartesian sensitivity set, not confidence coverage of a population response. Factor-two and the fixed fold grid are operational summaries, not authenticated biological significance limits. A positive lower bound excludes equality of the observed marker-state responses conditional on those envelopes. It does not isolate the intended edit or demonstrate protein preservation, full realised haplotype, isoform mechanism or external validation. Shared order is only a necessary condition for a common sequence-additive latent model, not a positive additive representation certificate. A feasible result is non-exclusion.

The final immutable freeze must include the current scientific module, runner, protocol and eligibility and all four declared source inputs. The launcher checks the entries supplied by that trusted root-owned freeze, not a cryptographic proof that a human granted authorization. No final real-data freeze/authorization exists among the inputs reviewed here. Outcome opening remains the root's separate gated action. A source-reproduction failure is recorded in JSON with no sensitivity summaries; downstream handling must inspect that status rather than infer scientific success from a normal process exit alone.

Source labels, three anonymous preparations and cross-cell comparisons do not establish independent cultures or paired states across markers/cells. Preserving within-row RPM pairing and avoiding invented cross-marker pairing is appropriate for this deliberately descriptive question. All missingness, finite-set scope and failed source gates must remain visible in the eventual result. Current corrected protocol supports these limits.

## Accepted revision and evidence

- Protocol: 5cb5711d220f91c179b8ee428dea0e91c72a725382032bc8a82083f0e5ec2e4f
- Core: 9d0bd0f3363e158e4df1cb028c7922b0bf94284bc6169138f894b072871c4cf9
- Runner: bfba04ed71fe695978327680e54fca45b197da42809db0632fe935a4e95d3cd7
- Preparation: 9ed008c98b8d51a6b7473404f12bb4e2321c326c666251d75021e816fa1e7cc9
- Eligibility: c5c2b20e6729918219c5cf3deb757ff242241cdf53f0f635f1dc7d677282ee32
- Author tests: 448b75c062da3fa52ef97a38629d4cd0f80a6c24cc9734e1baf19dfd5b9a68a8
- Author synthetic receipt: 90ac157ab715bb0eea70e8cd49650e402e45c372183466a726e2ba696d6d4d76

Independent scripts and receipts: verify_math.py / math_checks.json and verify_reader.py / reader_checks.json. They were separately authored instead of rerunning the author's test driver. This is an independent implementation challenge by the earlier interval-proof author, not an independent second derivation of that author's general theorem. No empirical source verification, empirical statistical calibration or new result is claimed.
