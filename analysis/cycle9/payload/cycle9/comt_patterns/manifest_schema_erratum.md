# Manifest schema correction

The first generic computation-manifest validation rejected the metadata adapter because it used output kind `evidence` rather than the allowed `log`, and omitted the declared result path. The original invalid manifest, initial final manifest and adapter source are preserved in `manifest_history_v1/`.

The corrected metadata uses `log` for the supporting evidence artifacts and declares `observed_v1/result.json`. That result filename was specified in the frozen `analyze.py` implementation before execution, and the output directory was supplied in the recorded invocation. This correction does not invent a new outcome declaration after inspecting results; it records the existing one in the validator's required field.

No scientific source, protocol, comparison, result, exact verifier or receipt changed. The once-only observed execution and distinct reconstruction were not rerun. This was a manifest-schema error, not an observed source-schema or arithmetic failure.
