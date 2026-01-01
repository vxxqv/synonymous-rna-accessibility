# Binary stored-float independent-reference correction

Status: post-result reference-only correction, ready for independent review. No new observed replay has run. Both the prior source-first reference and its failed exact comparison remain unchanged.

Two disclosed expectations change: finite Python floats are converted to their exact binary rational, matching the already-frozen production kernel; and the RNA-only linear model's irrelevant load_outside diagnostic is false, matching its constant-zero load input. Raw decimal strings and existing fractions remain exact. Every count-normalization, optimization, prediction, scoring, fold and scientific decision algorithm is unchanged. No biological quantity or logarithm is claimed intrinsically rational.

The complete code diff and change_inventory.json distinguish the two computational corrections from predecessor/amendment provenance and receipt-label updates. The first correction affects model inputs and observed target arithmetic through the common rational function, rather than rounding final answers or relaxing comparisons. The diagnostic correction changes no prediction.

## Invented evidence

- 46 direct non-dyadic fit cases compare920 predictions across all five models and92 coefficient/rank sets with the frozen production executable. These include rank-degenerate and nonnegative-slope boundary cases.
- A full invented125-variant/53-site CV comparison matches2,500 predictions,2,500 errors,80 aggregate scores,3,160 site errors,40 coefficient/rank sets and500 range-flag sets across both maps and both load types.
- Separate constant-response and explicit-missingness CV cases preserve zero-baseline failures and unavailable states, with every compared output exact.
- A load-extrapolation control has RNA-only load_outside=false and joint load_outside=true. Invalid RNA and load flags are still rejected. A1e-40 exact prediction perturbation is rejected despite identical float display.
- Raw decimal strings, existing fractions, invalid finite-number domains and32 metadata/provenance guards pass. -O/-OO are refused.

The first authorization-test attempt exceeded the Windows legacy path length in its invented temporary directory. No scientific call ran. The test now uses an extended-path temporary parent, checks containment before automatic cleanup, and passes; the initial directory is retained. This changes only the test harness, not the checker.

## Future authorization, after review

```text
python reference_binary_float.py --authorization EXPLICIT_AUTHORIZATION.json --out NEW_DIRECTORY
```

Use schema COMT_independent_binary_float_authorization_v3, status explicitly_authorized and the sealed reference_sha256. Supply the21 roles documented by the previous source-first reference, plus these four, making25 explicit observed-input roles:

| Additional role | Artifact |
|---|---|
| source_first_reference | Preserved reference_source_first.py, SHA222bffad... |
| source_first_preaccess_seal | Preserved successor_preaccess_seal.json, SHA8f9fc728... |
| original_replay_failure_seal | Preserved actual_v3_failure_seal.json, SHA3fcb9018... |
| reference_contract_amendment | This directory's reference_contract_amendment.md |

All predecessor, stage, amendment and selected payload hashes are required before numerical decoding. The full raw reconstruction cannot be downgraded to a downstream-only observed check. Original guards, exact comparisons, numerical-identity tolerance2e-12, invalid-denominator failure recording and no-retuning rule remain in force. No archive, source simulation, network, Git or installation is needed by the checker.

This amendment resolves reference expectations; it does not establish full observed agreement. Any newly authorized complete comparison must preserve discrepancies. Extraction/RF masks remain conditioned inputs, condition numbers are not independently reproduced, and stored hash lineage does not authenticate historical access order. Even full computational agreement would not establish physical protein amount, causal load information, kinetic compensation, calibrated uncertainty or assay equivalence.
