# Review of the Windows archive path fix

**Accepted for the bounded release-path repair.** Final reviewed `analysis/verify_bundle.py` SHA256 is `6d10f957cb1560b919613c10c990a027ef555b50f1cca9d158533e877c123386`.

The change adds `os` and normalizes the repository root to a Windows extended path. It does not omit any path or weaken a hash check. Independent AST comparison against commit `45ee5651dce6e668a17e39a956d8c1f673b41d09` confirms that `require`, `verify` and the complete main entrypoint are unchanged. All scientific checks, file enumeration, Git blob hashing, failure behavior and output handling remain the same.

The smoke check independently enumerated the full 5,640-file Git tree using both the ordinary and extended root as the working directory. The two inventories were byte-identical. Ordinary Windows paths falsely reported 32 files absent; the corrected root found every file. The longest tracked path is 292 characters. Its ordinary `is_file()` returned false, its extended-path check returned true, and its bytes reproduce the exact committed Git blob `a816893bc45ca764c9d1f50c08be43363ab5163b`.

Five prefix cases passed: drive-letter root, already-extended drive root, non-Windows guard, UNC root, and already-extended UNC root. The first candidate incorrectly prefixed UNC roots; this review identified that case and the root agent corrected it before the final smoke check. The final code uses `\\?\UNC\` plus the server/share portion for a UNC root and `\\?\` for a drive-letter root. No network share was contacted.

## Acceptance and falsifiers

Acceptance requires unchanged verification functions, the same full tracked inventory, access to every tracked path, a matching longest-file blob, and correct platform/prefix behavior. All these checks passed at the source hash above. A changed scientific AST, omitted path, mismatched blob, double prefix, malformed UNC root or altered non-Windows root would have failed this review.

This is a source and filesystem smoke check, not a rerun of the scientific verification. The complete required integrity check must still run against the committed final revision. A release archive extracted without Git needs the separate portable numerical replay; this existing exact-tree checker intentionally still requires Git.

Evidence: `check.py` contains the independent test, and `result_v1.json` records the final source hash, all five prefix cases, unchanged AST results, inventory count and longest-file blob comparison. No repository file was edited by this review.
