"""Numerical replay of the already reviewed PERSIST independent reconstruction.

Uses repository-relative archived inputs and an external fresh output directory.
No refolding, source-table parsing, Git, network or installation is performed.
"""
from __future__ import annotations

# This guard precedes imports, argument processing, dependency loading and all
# scientific input reads. The preserved synthetic checks use assertions.
if not __debug__:
    raise RuntimeError("Assertions must be enabled; -O and -OO are not permitted")

from datetime import datetime, timezone
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import platform
import subprocess
import sys
import traceback
import zipfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
EVALUATOR_SHA = "0711ae066da3c443312451e42fcbb651fad8280e5d7173838c3787b3990a1bd0"
PROTOCOL_SHA = "a3dc686a23acaa12af2f2f44154b8baac819cb6df3a0a6d81328dd3ebdc5fd24"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def relative(name):
    value = PurePosixPath(name)
    require(bool(value.parts) and not value.is_absolute() and "\\" not in name,
            f"Invalid archive-relative path: {name!r}")
    require(all(p not in (".", "..") and ":" not in p for p in value.parts), f"Unsafe path: {name!r}")
    return Path(*value.parts)


def file_at(repo, name):
    path = (repo / relative(name)).resolve()
    require(path.is_relative_to(repo) and path.is_file(), f"Missing or escaping repository input: {name}")
    return path


def input_hashes(repo, manifest):
    return {name: sha(file_at(repo, name)) for name in manifest["scientific_and_package_input_sha256"]}


def load_evaluator(path):
    spec = importlib.util.spec_from_file_location("preserved_persist_independent", path)
    require(spec is not None and spec.loader is not None, "Cannot load held evaluator")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def extract_evidence(source, destination, expected):
    destination.mkdir()
    with zipfile.ZipFile(source) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)) and set(names) == set(expected), "Historical evidence inventory mismatch")
        for info in archive.infolist():
            rel = relative(info.filename)
            require(not info.is_dir() and ((info.external_attr >> 16) & 0o170000) != 0o120000,
                    "Unsupported historical evidence member")
            path = destination / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(archive.read(info))
    observed = {name: sha(destination / relative(name)) for name in expected}
    require(observed == expected, "Historical evidence byte mismatch")
    return observed


def execute(out, repo, manifest, before):
    roles = {name: file_at(repo, path) for name, path in manifest["input_roles"].items()}
    preserved = out / "preserved_history"
    historical_hashes = extract_evidence(roles["historical_evidence_zip"], preserved, manifest["historical_evidence_sha256"])
    frozen = roles["evaluator"]
    require(sha(frozen) == EVALUATOR_SHA and sha(roles["protocol"]) == PROTOCOL_SHA, "Wrong frozen numerical source or protocol")
    # Preserve the old authorization as history. This local mapping reuses only
    # already authorized archived cells and never follows its old machine paths.
    authorization = {
        "status": "authorized_by_root", "operation": "independent_PERSIST_prediction_v2_compare",
        "independent_evaluator_sha256": EVALUATOR_SHA, "protocol_sha256": PROTOCOL_SHA,
        "artifacts": {name: {"path": str(roles[name]), "sha256": sha(roles[name])} for name in
                      ("features", "cds", "selected_targets", "author_predictions", "author_validation")},
        "portable_mapping_only": True,
        "authorization_semantics": "Relative remapping of historically authorized and archived 37 selected endpoint cells; not a new historical access-order claim",
        "historical_authorization_sha256": manifest["historical_authorization_sha256"]}
    mapping = out / "portable_authorization.json"
    write(mapping, authorization)
    jobs = [
        ("synthetic", [sys.executable, "-I", "-B", str(roles["synthetic_tests"]), "--protocol", str(roles["protocol"]), "--out", str(out / "synthetic")]),
        ("numerical", [sys.executable, "-I", "-B", str(frozen), "--protocol", str(roles["protocol"]),
                       "--authorization", str(mapping), "--out", str(out / "numerical")]),
    ]
    for name, command in jobs:
        completed = subprocess.run(command, cwd=out, text=True, encoding="utf-8", capture_output=True, timeout=240, check=False)
        (out / (name + "_stdout.txt")).write_text(completed.stdout, encoding="utf-8")
        (out / (name + "_stderr.txt")).write_text(completed.stderr, encoding="utf-8")
        require(completed.returncode == 0, f"{name} replay failed with exit {completed.returncode}; preserve this output")
    synthetic = read(out / "synthetic/receipt.json")
    numerical = read(out / "numerical/receipt.json")
    require(synthetic["status"] == "passed" and synthetic["tests_run"] == 12, "Frozen synthetic checks did not pass")
    require(numerical["status"] == "passed" and numerical["all_inputs_unchanged"] is True, "Frozen numerical comparison failed")
    evaluator = load_evaluator(frozen)
    targets = read(roles["selected_targets"])
    scale = max(1.0, *(abs(float(v)) for v in targets["development_targets"].values()),
                *(abs(float(v)) for v in targets["holdout_targets"].values()))
    atol = evaluator.ABSOLUTE_RESPONSE_SCALE_FACTOR * scale
    comparisons = {}
    raw_identity = {}
    for name in ("independent_predictions.json", "independent_validation.json"):
        expected = preserved / "independent_result_review/observed_v1" / name
        actual = out / "numerical" / name
        comparisons[name] = evaluator.compare_numeric(read(expected), read(actual), atol, name)
        raw_identity[name] = sha(expected) == sha(actual)
    original_synthetic = read(preserved / "independent_result_review/synthetic_v2/receipt.json")
    require(synthetic["counts"] == original_synthetic["counts"] and synthetic["tests_run"] == original_synthetic["tests_run"], "Historical synthetic coverage differs")
    current = read(out / "numerical/independent_validation.json")
    require(len(current["models"]) == 6 and len(current["rows"]) == 15
            and len(current["leave_one_holdout_row_out"]) == 15, "Incomplete fixed results")
    require(current["at_least_ten_percent_better_than_every_baseline"] is False,
            "Frozen failed primary decision changed")
    require(all(row["at_least_ten_percent_better_than_every_baseline"] is False
                for row in current["leave_one_holdout_row_out"]), "Frozen influence decision changed")
    after = input_hashes(repo, manifest)
    require(before == after, "Repository/archive input changed during replay")
    require(not list(out.rglob("__pycache__")), "Unexpected generated bytecode")
    generated = ["portable_authorization.json", "synthetic/receipt.json", "numerical/receipt.json",
                 "numerical/independent_predictions.json", "numerical/independent_validation.json"]
    return {
        "status": "passed", "all_six_models_and_fifteen_influence_records_reproduced": True,
        "synthetic_tests": 12, "historical_output_comparisons": comparisons,
        "raw_independent_output_bytes_identical_on_this_host": raw_identity,
        "author_comparison_receipt": "numerical/receipt.json", "primary_success": False,
        "all_fifteen_influence_success_decisions": False,
        "input_sha256_before": before, "input_sha256_after": after, "all_inputs_immutable": True,
        "historical_evidence_sha256": historical_hashes,
        "generated_output_sha256": {name: sha(out / relative(name)) for name in generated},
        "historical_git_or_access_order_authenticated": False, "refolding_performed": False,
        "combined_source_table_read": False, "runtime_network_or_installation_used": False,
        "scope": "Numerical replay from fixed MUP/features and 37 selected endpoint cells; no cross-platform folding, causal or calibrated-uncertainty claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path, help="Extracted computational repository or source archive root")
    parser.add_argument("--out", required=True, type=Path, help="New directory outside the repository")
    args = parser.parse_args()
    repo, out = args.repo.resolve(), args.out.resolve()
    require(repo.is_dir(), "Repository/source archive root does not exist")
    require(not out.exists(), "Refusing to overwrite existing output")
    require(not out.is_relative_to(repo), "Output must be outside the repository")
    manifest = read(HERE / "hash_manifest.json")
    before = input_hashes(repo, manifest)
    require(before == manifest["scientific_and_package_input_sha256"], "Frozen package/reference input hash mismatch")
    import numpy as np
    require(np.__version__ == "2.5.4", "Preinstalled NumPy 2.5.4 is required; no automatic installation")
    out.mkdir(parents=True)
    receipt = {"schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
               "software": {"Python": platform.python_version(), "NumPy": np.__version__},
               "hardware": platform.platform() + "; " + platform.machine(),
               "wrapper_sha256": sha(Path(__file__)), "input_manifest_sha256": sha(HERE / "hash_manifest.json"),
               "output_outside_repository": True, "assertions_enabled": __debug__}
    try:
        receipt.update(execute(out, repo, manifest, before))
    except Exception as error:
        receipt.update({"status": "failed", "error": f"{type(error).__name__}: {error}", "retuned": False})
        (out / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write(out / "receipt.json", receipt)
        print(json.dumps({"status": "failed", "error": receipt["error"]}))
        return 1
    write(out / "receipt.json", receipt)
    print(json.dumps({"status": "passed", "models": 6, "influence_records": 15,
                      "primary_success": False, "historical_access_order_authenticated": False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
