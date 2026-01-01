# Location of preserved originals

The original `ERRATUM.md` is preserved as supplied by the independent reviewer. Its reference to `history_v1` describes the reviewer's source workspace. In this computational checkout, the original verifier, result and receipt are retained unchanged in the parent `supplement_review` directory. The corrected files are in `correction_v2`. No separate `history_v1` directory is required here.

The original and corrected source hashes, result hashes and the scope-only comparison are pinned in `erratum_receipt.json`. Portable execution uses `analysis/cycle5/replay_supplement_reviews.py`; the historical `run_correction.py` retains the review workspace layout.
