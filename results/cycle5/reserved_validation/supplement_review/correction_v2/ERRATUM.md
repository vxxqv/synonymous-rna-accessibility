# Scope wording correction

The first reserved review carried a copied sentence describing development-only access and denying reserved-outcome access. That sentence was incorrect. The review actually used the explicitly authorized reserved evaluation and diagnostics, as its arm, input hashes, calculations and original receipt already recorded.

Only the scope string was changed in the copied verifier. All original package files are retained unchanged, with byte-identical copies under history_v1. The corrected verifier was rerun on identical pinned inputs. Only scope and the self-reported reviewer-source hash differ in the result; every scientific field and check count is exactly identical. All 13 numerical/helper function ASTs are unchanged. The separate runner adjusts workspace discovery for the extra package depth; it does not alter inputs, models, thresholds or calculations.

The correction does not change the PASS arithmetic verdict, failed coverage or failed prediction screen.
