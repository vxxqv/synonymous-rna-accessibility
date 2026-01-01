# PERSIST independent numerical replay

This directory replays the already frozen independent numerical reconstruction. It does not open the full source table, refold RNA, fit an alternative model, retune a grid or select new outcomes. All six fixed models, ridge tuning grids and OOF predictions, primary/secondary metrics and all 15 influence records are checked. The failed primary result is preserved.

## Requirements and command

Use Python 3.12 with NumPy **2.5.4 preinstalled**. The verified run used CPython 3.12.14. `requirements.txt` records the sole external dependency. No installation, Git query or network call occurs at runtime. ViennaRNA, Biopython, NetworkX and the original Windows layout are not required.

From the repository or extracted computational source archive, run:

```text
python -I -B analysis/programme/retrospective/persist_v2/independent_replay/replay.py --repo . --out ../persist-independent-replay
```

The output path must be new and outside the repository. An existing path or an in-repository destination is rejected. `-O` and `-OO` are rejected before argument processing, scientific input reads, dependency imports or output creation. The wrapper and frozen sources also disable bytecode writing. Run in a source archive without `.git` if desired. A normal run takes seconds; each subprocess has a 240-second timeout. Failures are preserved and return a nonzero exit status. No mathematical source is repaired or altered by replay.

## Referenced and included inputs

Existing archived inputs are referenced without duplication:

- `analysis/programme/retrospective/persist_v2/protocol.json` and the author interface source `persist_prediction.py`.
- `results/programme/persist_v2/train/features.json` and `predictions.json`.
- `results/programme/persist_v2/validate/validation.json`.

New input subsets are under `results/programme/persist_v2/independent_replay/inputs/`: `selected_targets.json` contains exactly the authorized endpoint for all 22 development and 15 holdout rows; `cds.json` contains only CDS metadata for those 37 rows. No historical AUP or other outcome column is included. Both subset files are byte-preserved from the root-authorized numerical audit.

`frozen/evaluator.py` is the exact pre-outcome independent source, SHA256 `0711ae066da3c443312451e42fcbb651fad8280e5d7173838c3787b3990a1bd0`. Its unchanged synthetic checker is beside it. `results/programme/persist_v2/independent_replay/original_evidence.zip` preserves the original source, freeze, history, design, exposure ledger, synthetic receipts, observed independent predictions/evaluation, review and historical authorization. Their individual original byte hashes are in `hash_manifest.json`. Old absolute paths are inactive historical provenance, never runtime dependencies.

`hash_manifest.json` pins every scientific input and static package component. `integration_manifest.json` lists all integrated files, including the successful portable verification receipt, and identifies its own usual self-hash exclusion. The release ZIP's externally supplied SHA256 anchors that complete manifest.

## What runs and what is compared

The wrapper verifies every input hash before decoding numeric results. It safely extracts the historical evidence into the external output folder. It creates an explicit local path mapping for the five already authorized numerical artifacts, with the original evaluator/protocol pins. The mapping's `authorized_by_root` marker is the frozen interface's required historical contract value; the additional metadata states that it is only a portable remapping, not a newly authenticated access event.

The unchanged 12 synthetic checks run first. Then the unchanged evaluator reconstructs all six predictions using augmented SVD ridge, the fixed 15 development groups, training-only transforms and the original tie rule. It writes independent answers before decoding author numeric answers for comparison. The fresh complete independent prediction and evaluation payloads are also compared against the preserved independent audit. No frozen scientific function is edited or wrapped with a different mathematical method.

The original predeclared comparison tolerances apply: response-unit `atol=1e-10*max(1,max|y_dev|,max|y_holdout|)`, `rtol=1e-9`; dimensionless rank and gain metrics use `atol=1e-12, rtol=1e-10`. Selected penalties, nearest-CDS ties/distances, identifier order, success decisions and undefined statuses must agree exactly. The fresh evaluator's original receipt details all author comparisons. The portable receipt also reports raw independent output byte identity on the tested host. Bitwise identity across BLAS libraries or platforms is not claimed.

Successful execution writes `receipt.json`, both independent JSON payloads under `numerical/`, the original frozen evaluator's receipt, input hashes before/after, original evidence hashes, and synthetic check evidence. Every output is outside the repository. Historical runtime, absolute commands, software metadata and old Git fields are retained as history rather than being recertified. Source archives need no historical Git objects to run the numerical comparison.

## Result and limits

The original failed primary decision remains false. Primary MAE is approximately 0.031442, versus 0.030139 for composition alone and 0.031479 for the codon/composition control. Every one-row omission also fails the fixed complete-baseline rule. These finite comparisons do not justify outcome-driven model or cohort revision.

MUP and all feature columns are fixed archived numerical inputs. This replay does not independently authenticate folding binaries, reproduce cross-platform partition functions, recover the raw target table, establish historical access order or verify Git chronology. The selected-target extraction and prior access controls remain recorded provenance. It provides numerical reproduction on one selected reporter/context/study, without a confidence interval, causal mechanism, calibrated rate or half-life, gene-level independence, or generalization to another gene. No PDF, manuscript, figure, unrelated outcome, dependency tree or bytecode cache belongs to this package.
