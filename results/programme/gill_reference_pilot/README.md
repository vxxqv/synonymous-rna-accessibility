# Frozen reference-protein RNA pilot

The predeclared predictive criterion failed. At pseudocount 0.5, pairwise shrinkage has held-out RMSE 0.5506172851, versus 0.4656950429 for the constant comparator: 18.2356% higher error. All seven predictors, all 32 genotypes and all three pseudocounts remain in the result tables. The target also fails at 0.25 and 1. No model or endpoint is retuned on these examined outcomes.

The conclusion concerns the specified procedure on one held-out RNA library, conditional on the already inspected DNA denominator. The complete K2 x K4 x K4 design makes numerical interaction contrasts identifiable, but unknown biological preparation relationships and measurement error prevent a calibrated biological-interaction conclusion. Distinct deposited files and sample identifiers do not establish independent cultures. Nine RNA-by-DNA contrasts are dependent sensitivity scenarios, not nine biological validation samples.

## Frozen chronology

- Protocol and reviewed implementation: commit `faa9bb248ef3969495f54504987f63d46dcad133`, before any RNA fitting.
- All model predictions and sensitivities: commit `d186707a74a40e14a3ed9175337796d9de25c3dd`, before CDS3 RNA was decoded.
- Prediction SHA256: `ff04894c283310fc9f47ae01f35191bc6a9be83493e39baaaf29d6da9bdd3660`.
- `prevalidation_gate.json` independently compares the exact committed prediction, protocol and source bytes with local bytes and verifies the remote branch before validation.
- The validation RNA library is now examined and must not be relabelled untouched or used for confirmatory retuning. Outcomes from other protein fibres were not decoded.

## Reproduction

Use the pinned repository environment. Retrieve the original processed CDS workbook from [Mendeley Data version 1](https://doi.org/10.17632/ds7zymk875.1). Its required SHA256 is `0f80ba9868c54548db938c65d09d8caa0419c5721bfa379949f6a18ad2568e89`; the workflow rejects any other bytes. Source provenance is in `results/programme/retrospective/` and the six-library mapping in `analysis/programme/retrospective/gill_replicate_metadata.json`. The deposited parser and construct audits establish the intended mapping; raw read processing has not been reproduced here.

```text
python analysis/programme/retrospective/gill_reference_pilot.py synthetic --outdir work/gill_replay/synthetic
python analysis/programme/retrospective/gill_reference_pilot.py train --protocol analysis/programme/retrospective/gill_reference_protocol_v1.json --eligibility results/programme/retrospective/gill_amplicon_pilot_v2.json --workbook /path/to/cds_processed.xlsx --outdir work/gill_replay/training
python analysis/programme/retrospective/gill_reference_pilot.py validate --protocol analysis/programme/retrospective/gill_reference_protocol_v1.json --eligibility results/programme/retrospective/gill_amplicon_pilot_v2.json --workbook /path/to/cds_processed.xlsx --predictions results/programme/gill_reference_pilot/training/predictions.json --prediction-sha256 ff04894c283310fc9f47ae01f35191bc6a9be83493e39baaaf29d6da9bdd3660 --outdir work/gill_replay/validation
python -m pytest tests/test_gill_reference_pilot.py -q
```

Output directories must be new. The validation command uses the historically frozen predictions, so reproduction does not silently create a replacement prediction freeze. Recomputed training predictions may be checked against them. Environment, runtime and repository-status fields can differ on replay; scientific predictions and losses must agree. Reproduction after outcome exposure is not new prospective validation.

The `training` and `validation` directories include complete selected count vectors, predictions, all errors, all 496 endpoint pairs, critical-threshold summaries and labelled sensitivity scenarios. These small derived tables permit independent arithmetic replay without downloading raw sequencing reads. The reference anchor is fixed. Centered RMSE is only a predeclared reference-error diagnostic; no confidence intervals or p-values are provided. Positive observed-field barriers do not establish a noise-free biological barrier.

### Independent saved-vector replay

The default independent result audit requires a full Git clone containing both frozen commits. It verifies their ancestry and committed blobs as well as every stored input/output hash and numerical result:

```text
python analysis/programme/retrospective/gill_result_independent.py --outdir work/gill-independent-full-history
```

For a source archive without `.git`, explicitly request arithmetic and stored-hash replay only:

```text
python analysis/programme/retrospective/gill_result_independent.py --source-only --outdir work/gill-independent-source-only
```

Source-only mode makes no Git query, records `historical_git_verified=false`, and cannot verify historical commit or access-order evidence. It retains every arithmetic and stored-hash check and never opens the observed workbook. Its conservative dirty-tree sentinel means “not certified clean,” not an assessed Git state. Both modes require a fresh output directory inside repository `work/` or `results/`. The default fails when the required history is absent; it never silently falls back.

The amended auditor's full-history receipt is `independent_audit/portable/`; `independent_audit/current/` and the earlier receipts are historical records of earlier auditor versions. The pre-portability source is preserved under `independent_audit/history/`. The portable receipt also records the no-history archive-mode checks and exact equality of scientific results between modes.

## Source attribution and terms

Manraj Gill, *Data for: Long-range mRNA folding shapes expression and sequence of bacterial genes*, Mendeley Data, version 1 (2026), DOI 10.17632/ds7zymk875.1, CC BY 4.0. The associated paper is Gill, Kim, Xue, Thappeta, Taggart and Li, *Molecular Cell* 86(18), 3789 to 3802.e11 (2026), DOI 10.1016/j.molcel.2026.08.014. The source already establishes its biological folding mechanism; this pilot does not claim that discovery.

The retained counts are an attributed subset of the source deposit. This reanalysis changes the response estimator, normalization and prediction models. Source data retain [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and are not silently relicensed by the repository's software license. No source-author endorsement is implied. The computational protocol, selective reader, synthetic controls and result summaries are provided for transparent replication of the failed comparison.
