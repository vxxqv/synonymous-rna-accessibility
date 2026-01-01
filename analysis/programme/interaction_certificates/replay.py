"""Replay the byte-preserved A3 finite checks without Git or the original workspace.

Only filesystem/provenance orchestration is adapted. Mathematical functions are
imported from the frozen payload without source rewriting.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import platform
import shutil
import subprocess
import sys
from time import perf_counter
import traceback
import zipfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PRODUCT_SHA256 = "96c6d35a6243dcf1ba0dbd08100d82dbbf553a583e3844e9ac6a31646d5ec63a"
AUTHOR_CORE = ("named_examples", "finite_suite", "maxcut_reduction_checks", "counts")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def safe_relative(name: str) -> Path:
    # Historical freeze keys use Windows separators. Interpret them portably;
    # do not rewrite the preserved JSON bytes or use historical absolute paths.
    pure = PurePosixPath(name.replace("\\", "/"))
    require(bool(pure.parts) and not pure.is_absolute(), f"Invalid relative path: {name!r}")
    require(all(p not in (".", "..") and ":" not in p for p in pure.parts), f"Unsafe path: {name!r}")
    return Path(*pure.parts)


def snapshot(root: Path, names) -> dict[str, str]:
    return {name: sha(root / safe_relative(name)) for name in names}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"Cannot load {path.name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def canonical(value):
    return json.loads(json.dumps(value, default=str))


def first_difference(expected, actual, at="$") -> str | None:
    if type(expected) is not type(actual):
        return f"{at}: {type(expected).__name__} != {type(actual).__name__}"
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            return f"{at}: different key sets"
        for key in expected:
            difference = first_difference(expected[key], actual[key], f"{at}.{key}")
            if difference:
                return difference
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            return f"{at}: different list lengths"
        for i, (left, right) in enumerate(zip(expected, actual)):
            difference = first_difference(left, right, f"{at}[{i}]")
            if difference:
                return difference
    elif expected != actual:
        return f"{at}: {expected!r} != {actual!r}"
    return None


def compare(expected, actual, ignored=()) -> dict:
    # Only explicitly named top-level timing/provenance fields are excluded.
    # Counters, statuses, exact values, qualifications and list order all count.
    require(set(ignored) <= set(expected), "Unknown ignored field in historical receipt")
    require(set(expected) == set(actual), "Different receipt top-level key sets")
    left = {k: v for k, v in expected.items() if k not in ignored}
    right = {k: v for k, v in actual.items() if k not in ignored}
    difference = first_difference(left, right)
    require(difference is None, f"Exact result comparison failed: {difference}")
    return {"exact_agreement": True, "ignored_top_level_fields": list(ignored),
            "compared_top_level_fields": sorted(left)}


def validate_payload(manifest: dict, preserved: Path) -> dict:
    names = manifest["payload_sha256"]
    observed = snapshot(preserved, names)
    require(observed == names, "Frozen payload byte hashes differ")
    product = "outputs/synonymous-rna-accessibility/src/product_closure.py"
    require(observed[product] == PRODUCT_SHA256, "Wrong additive dependency")
    freeze_checks = {}
    for rel in manifest["freeze_files"]:
        frozen = read(preserved / safe_relative(rel))
        inputs = frozen["input_sha256"]
        for name, expected in inputs.items():
            path = safe_relative(name).as_posix()
            require(path in names, f"Frozen dependency absent from payload: {path}")
            require(observed[path] == expected, f"Historical freeze mismatch: {path}")
        freeze_checks[rel] = {"verified_input_count": len(inputs), "exact_hash_match": True}
    return freeze_checks


def replay(out: Path, manifest: dict) -> dict:
    preserved, reconstruction = out / "preserved", out / "reconstruction"
    preserved.mkdir()
    reconstruction.mkdir()
    with zipfile.ZipFile(HERE / "computational_sources.zip") as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), "Duplicate ZIP members")
        require(set(names) == set(manifest["payload_sha256"]), "ZIP contents differ from manifest")
        for info in archive.infolist():
            rel = safe_relative(info.filename)
            require(not info.is_dir() and ((info.external_attr >> 16) & 0o170000) != 0o120000,
                    "Directories and symlinks are not payload files")
            destination = preserved / rel
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(archive.read(info))
    freeze_checks = validate_payload(manifest, preserved)
    before = snapshot(preserved, manifest["payload_sha256"])
    # Generated historical results/manifests stay in preserved/. Only static
    # inputs go to the active layout, so the original reference cannot overwrite
    # an archived receipt and the creator's run starts without results.json.
    for name in manifest["active_input_paths"]:
        rel = safe_relative(name)
        destination = reconstruction / rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(preserved / rel, destination)
    active_before = snapshot(reconstruction, manifest["active_input_paths"])
    base = reconstruction / "work/interaction_cycle3"
    expected_base = preserved / "work/interaction_cycle3"
    frozen = read(base / "freeze.json")
    require(not (base / "results.json").exists(), "Author output already exists")

    started = perf_counter()
    with (out / "author_stdout.txt").open("w", encoding="utf-8") as stdout, \
            (out / "author_stderr.txt").open("w", encoding="utf-8") as stderr, \
            redirect_stdout(stdout), redirect_stderr(stderr):
        author = load_module("a3_portable_author", base / "verify_interactions.py")
        result = {"named_examples": author.named_examples(),
                  "finite_suite": author.verify_suite(),
                  "maxcut_reduction_checks": author.maxcut_checks(),
                  "counts": dict(author.COUNTS)}
    author_seconds = perf_counter() - started
    result.update({"runtime_seconds": author_seconds, "software": {"Python": platform.python_version()},
                   "status": "passed_for_frozen_finite_inputs",
                   "proof_scope": "Finite implementation corroboration; universal assertions rely on claim_and_proof.md",
                   "input_sha256": frozen["input_sha256"]})
    result = canonical(result)
    write(base / "results.json", result)
    author_comparison = compare(read(expected_base / "results.json"), result,
                                ("runtime_seconds", "software"))
    legacy = read(expected_base / "replay_results.json")
    require(set(legacy) == set(AUTHOR_CORE), "Unexpected historical replay receipt fields")
    author_comparison["prior_full_replay_core_agreement"] = compare(
        legacy, {k: result[k] for k in AUTHOR_CORE})["exact_agreement"]

    independent = base / "independent_review"
    (independent / "run").mkdir()
    started = perf_counter()
    # The reference has its own implementation and executes as a separate
    # process. Its only creator imports are clearly labelled storage diagnostics.
    completed = subprocess.run([sys.executable, "-B", str(independent / "reference.py")],
                               cwd=reconstruction, capture_output=True, text=True,
                               encoding="utf-8", timeout=240, check=False)
    independent_seconds = perf_counter() - started
    (out / "reference_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (out / "reference_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    require(completed.returncode == 0, f"Independent reference exited {completed.returncode}")
    reference_result = read(independent / "run/results.json")
    reference_comparison = compare(read(expected_base / "independent_review/run/results.json"), reference_result)

    started = perf_counter()
    with (out / "revision2_stdout.txt").open("w", encoding="utf-8") as stdout, \
            (out / "revision2_stderr.txt").open("w", encoding="utf-8") as stderr, \
            redirect_stdout(stdout), redirect_stderr(stderr):
        corrected = load_module("a3_portable_corrected", base / "revision_2/targeted_regression_v2.py")
        require(sha(corrected.SOURCE) == corrected.CODE_SHA, "Correction imports wrong author source")
        require(sha(corrected.REVIEW) == corrected.REVIEW_SHA, "Correction cites wrong independent review")
        fresh_author = load_module("a3_portable_accounting_author", base / "verify_interactions.py")
        corrected_result = corrected.suite(fresh_author)
    corrected_seconds = perf_counter() - started
    corrected_result.update({"status": "targeted_regressions_passed", "runtime_seconds": corrected_seconds,
        "provenance": "Portable call to unchanged revision-2 suite; historical Git queries are not repeated."})
    corrected_result = canonical(corrected_result)
    write(base / "revision_2/results_v2.json", corrected_result)
    correction_comparison = compare(read(expected_base / "revision_2/results_v2.json"), corrected_result,
                                    ("runtime_seconds", "provenance"))

    after = snapshot(preserved, manifest["payload_sha256"])
    active_after = snapshot(reconstruction, manifest["active_input_paths"])
    require(before == after and active_before == active_after, "An input changed during replay")
    require(not list(out.rglob("__pycache__")), "Unexpected bytecode cache in replay")
    require(not list(out.rglob(".git")), "Unexpected Git metadata in replay")
    result_paths = ["work/interaction_cycle3/results.json",
                    "work/interaction_cycle3/independent_review/run/results.json",
                    "work/interaction_cycle3/revision_2/results_v2.json"]
    return {"status": "passed", "scientific_agreement": {
                "author_finite_receipt": author_comparison,
                "separately_authored_reference": reference_comparison,
                "corrected_revision2_diagnostics": correction_comparison},
            "historical_git_verified": False,
            "historical_git_note": "Original receipts are preserved byte-for-byte; their Git provenance was not reauthenticated.",
            "no_git_metadata_in_replay": True, "frozen_inputs_unchanged": True,
            "payload_sha256_before": before, "payload_sha256_after": after,
            "active_input_sha256_before": active_before, "active_input_sha256_after": active_after,
            "historical_freeze_checks": freeze_checks,
            "generated_result_sha256": snapshot(reconstruction, result_paths),
            "runtime_seconds": {"author": author_seconds, "independent_reference": independent_seconds,
                                "revision2": corrected_seconds},
            "finite_coverage": {"author_parameter_inputs": 2818,
                 "author_maxcut_graphs": len(result["maxcut_reduction_checks"]),
                 "reference_Q2_parameter_inputs": reference_result["Q2_fields"],
                 "reference_Q3_parameter_inputs": reference_result["Q3_rational_fields"],
                 "reference_projection_checks": reference_result["finite_projection_checks"],
                 "reference_scaled_maxcut_checks": len(reference_result["maxcut_checks"]),
                 "distinct_field_qualification": corrected_result["coverage_qualification"]},
            "current_contract": {
                "proof": "preserved/work/interaction_cycle3/revision_2/claim_and_proof.md",
                "complexity": "preserved/work/interaction_cycle3/revision_2/complexity_contract.md",
                "correction_ledger": "preserved/work/interaction_cycle3/revision_2/correction_ledger.md",
                "trace_erratum": "preserved/work/interaction_cycle3/revision_2/trace_erratum.md",
                "scope": "Finite corroboration only; no new generic method, calibrated uncertainty or biological result."}}


def main() -> int:
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not run with -O or -OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="New output directory; must not already exist")
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError("Refusing to overwrite an existing replay directory")
    manifest = read(HERE / "integration_manifest.json")
    for name, expected in manifest["package_sha256"].items():
        require(sha(HERE / safe_relative(name)) == expected, f"Package input differs: {name}")
    out.mkdir(parents=True)
    receipt = {"schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
               "software": {"Python": platform.python_version(), "implementation": platform.python_implementation()},
               "package_sha256": manifest["package_sha256"],
               "adapter_sha256": sha(Path(__file__)),
               "assertions_enabled": __debug__, "network_used": False, "git_required": False}
    try:
        receipt.update(replay(out, manifest))
    except Exception as exc:
        receipt.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        (out / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write(out / "receipt.json", receipt)
        print(json.dumps({"status": "failed", "error": receipt["error"]}))
        return 1
    write(out / "receipt.json", receipt)
    print(json.dumps({"status": "passed", "all_three_exact_comparisons": True,
                      "historical_git_verified": False, "receipt": "receipt.json"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
