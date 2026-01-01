# Independent headerless-format adapter review

**PASS for the inspected adapter.** SHA256 `11addfe6ee9c02f7fc17b9b6a414167cc0d98baccbab71bf7fa899fd7b713bb6`; amendment SHA256 `73a89257046695a1882c01b98971018d78aba95a144afbf5e6752a9e10d2b9a62`.

The21 synthetic checks show that the first native-reference row, zero-valued rows, large exact integer counts and every other valid row are retained. LF, CRLF, absent final newline, blank internal lines and ordinary whitespace give the same dictionaries as the unchanged frozen parser with the expected header. An instrumented parser call confirms that the adapter passes exactly `HEADER + original_bytes` and no edited data row. It does not run or reproduce an author importer that discards the first row.

The original Git blob is checked against the original bytes before decoding or invoking the parser. A wrong hash on invalid UTF-8, and a one-count synthetic edit with the original hash, both fail at the source check. Duplicate first or later genotypes, incorrect ancestor mutation annotations, negative/fractional/non-ASCII counts, malformed rows, outside-design sequences, empty/all-zero libraries, an unexpected existing header and a different first genotype are refused. A BOM-prefixed headerless file is also refused because the BOM would occur inside the prefixed stream. This is fail-closed; no broader encoding support is claimed.

All nine files listed by the original scientific freeze match both their SHA256 entries and their exact Git bytes at commit `7d8faf59fddde21ef4f1177f13ab5d2144d192e6`. The original freeze itself is byte-identical to that commit. No endpoint, normalization formula, graph, model, domain or scientific exclusion changed. The frozen source contract supplies twelve nonempty40-hex-character Git blob identifiers to the CLI parser. The optional unpinned parser mode is used only for synthetic comparisons in this audit.

At inspection, `format_freeze.json` was not yet present. Therefore this verdict covers the adapter source, original freeze preservation and synthetic parsing, not a completed production CLI run or final amendment freeze. Root must create and pin the amendment manifest; `main` then checks its listed file hashes before replacing the parser and calls the unchanged core CLI, which enforces the original scientific freeze and no-overwrite guard. No deposited count file or biological outcome was opened here.

The amendment correctly limits the secondary point calculation to **corrected headerless source-formula reconstruction**. Its retained legacy JSON label does not authenticate historical published point estimates. The reported released R import discrepancy does not establish that the authors used these same headerless bytes in their actual historical working pipeline.

Evidence: `review_adapter.py` and `receipt.json`, SHA256 `373892aa747161d267b7ce6f1ecacf2dfef07599f6b9aa6712af23ea000eec42`. The test is entirely synthetic and separately authored; no shared source was edited.
