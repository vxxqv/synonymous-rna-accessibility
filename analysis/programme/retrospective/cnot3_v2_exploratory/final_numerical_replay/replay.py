"""Offline archive replay wrapper. Standard library only; numerical verifier unchanged."""
import sys
if not __debug__:
    raise RuntimeError("Optimized Python is refused before input access")
sys.dont_write_bytecode = True

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess

MANIFEST_SHA256 = "379b1cdbe4fd394fa6137e45e2fe3eec9b3a0352cda3108d7c0d733ee523a2ff"
SCHEMA = "CNOT3_final_numerical_archive_replay_v1"
VERIFIER_ROLES = (
    "training_features", "training_targets", "holdout_features", "design_receipt",
    "predictions", "protocol", "author_code", "evaluation", "holdout_targets",
    "commit_receipt", "authorization",
)


def native(path):
    """Use explicit long paths on Windows, including UNC paths."""
    value = os.path.abspath(path)
    if os.name == "nt" and not value.startswith("\\\\?\\"):
        value = "\\\\?\\UNC\\" + value[2:] if value.startswith("\\\\") else "\\\\?\\" + value
    return Path(value)


def digest(path):
    h = hashlib.sha256()
    with native(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def strict_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate JSON key: " + key)
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError("Nonfinite JSON constant: " + value)
    return json.loads(native(path).read_text(encoding="utf-8"),
                      object_pairs_hook=unique, parse_constant=nonfinite)


def contained_file(root, relative):
    if not isinstance(relative, str) or "\\" in relative or ":" in relative:
        raise ValueError("Invalid archive relative path")
    parts = PurePosixPath(relative)
    if parts.is_absolute() or not parts.parts or any(x in (".", "..") for x in parts.parts):
        raise ValueError("Invalid archive relative path")
    resolved_root = native(root).resolve(strict=True)
    path = (resolved_root / parts).resolve(strict=True)
    if not path.is_relative_to(resolved_root) or not path.is_file():
        raise ValueError("Archive path escapes root or is not a file: " + relative)
    return path


def checked_inputs(root, manifest):
    paths, hashes = {}, {}
    for section in ("inputs", "supporting_evidence"):
        for role, item in manifest[section].items():
            key = section + "/" + role
            path = contained_file(root, item["path"])
            actual = digest(path)
            if actual != item["sha256"]:
                raise ValueError("SHA256 mismatch before scientific input decoding: " + key)
            hashes[key] = actual
            if section == "inputs":
                paths[role] = path
    return paths, hashes


def write_json(path, value):
    native(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def replay(archive_root, output):
    manifest_path = native(__file__).parent / "manifest.json"
    if digest(manifest_path) != MANIFEST_SHA256:
        raise ValueError("Replay manifest SHA256 mismatch")
    manifest = strict_json(manifest_path)
    if manifest["schema"] != SCHEMA:
        raise ValueError("Unsupported replay manifest schema")
    root = native(archive_root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("Archive root is not a directory")
    out = native(output).resolve(strict=False)
    if out.is_relative_to(root):
        raise ValueError("Replay output must be outside the archive")
    if out.exists():
        raise FileExistsError("Replay output must be a new directory")
    paths, hashes = checked_inputs(root, manifest)
    expected = strict_json(paths["expected_verification"])
    if expected["status"] != "passed" or expected["failures"]:
        raise ValueError("Archived independent verification did not pass")
    for key, value in manifest["expected"].items():
        if expected[key] != value:
            raise ValueError("Archived verification counter mismatch: " + key)
    if expected["runtime"] != manifest["runtime"]:
        raise ValueError("Archived verification runtime mismatch")
    out.mkdir(parents=True, exist_ok=False)
    try:
        stage = out / "unchanged_verifier"
        stage.mkdir()
        for role, filename in (("verifier", "verify.py"), ("reference_core", "reference_core.py")):
            shutil.copyfile(paths[role], stage / filename)
            if digest(stage / filename) != manifest["inputs"][role]["sha256"]:
                raise ValueError("Copied verifier SHA256 mismatch: " + role)
        command = [str(native(sys.executable)), "-I", "-B", str(stage / "verify.py")]
        for role in VERIFIER_ROLES:
            command.extend(["--" + role.replace("_", "-"), str(paths[role])])
        command.extend(["--out", str(out / "verification")])
        # List arguments avoid a shell. No Git, network, environment installation,
        # author imports, or candidate selection occurs in this wrapper.
        process = subprocess.run(command, cwd=out, capture_output=True, text=True, encoding="utf-8")
        (out / "stdout.txt").write_text(process.stdout, encoding="utf-8")
        (out / "stderr.txt").write_text(process.stderr, encoding="utf-8")
        if process.returncode:
            raise RuntimeError("Unchanged verifier refused or failed; inspect stderr.txt and verification/failure.json")
        result_path = out / "verification" / "verification.json"
        result = strict_json(result_path)
        # Reuse the verifier's published tolerances. The historical maximum
        # discrepancy is diagnostic, not a newly tightened acceptance threshold.
        required = ("status", "authorization_status", "input_sha256", "verifier_sha256",
                    "reference_core_sha256", "runtime", "failures", "tolerances", "limits")
        for key in required + tuple(manifest["expected"]):
            if result[key] != expected[key]:
                raise ValueError("Replay receipt disagreement: " + key)
        _, after = checked_inputs(root, manifest)
        if after != hashes:
            raise ValueError("An archive input changed during replay")
        for role, filename in (("verifier", "verify.py"), ("reference_core", "reference_core.py")):
            if digest(stage / filename) != manifest["inputs"][role]["sha256"]:
                raise ValueError("Copied verifier changed during replay: " + role)
        receipt = {
            "status": "passed", "schema": SCHEMA,
            "wrapper_sha256": digest(__file__), "manifest_sha256": MANIFEST_SHA256,
            "checked_archive_files": len(hashes), "archive_sha256": hashes,
            "all_archive_inputs_unchanged": True,
            "verifier_copied_byte_identically": True,
            "verification_sha256": digest(result_path),
            "archived_verification_sha256": manifest["inputs"]["expected_verification"]["sha256"],
            "exact_receipt_content_match": result == expected,
            "exact_receipt_bytes_match": digest(result_path) == digest(paths["expected_verification"]),
            "checks": result["checks"], "numerical_values_compared": result["numerical_values_compared"],
            "maximum_abs_difference": result["maximum_abs_difference"],
            "archived_maximum_abs_difference": expected["maximum_abs_difference"],
            "runtime": result["runtime"], "tolerances": result["tolerances"],
            "limits": expected["limits"],
            "scope": "Offline numerical replay only. Cohort review evidence hash-checked, not rerun. No historical Git verification or new experimental validation.",
        }
        write_json(out / "replay_receipt.json", receipt)
        return receipt
    except Exception as error:
        write_json(out / "replay_failure.json", {"status": "failed", "error": str(error),
                                                "type": type(error).__name__})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", required=True)
    parser.add_argument("--out", required=True, help="New directory outside the archive")
    args = parser.parse_args()
    result = replay(args.archive_root, args.out)
    print(json.dumps({k: result[k] for k in ("status", "checks", "numerical_values_compared",
                                           "maximum_abs_difference", "exact_receipt_bytes_match")}))


if __name__ == "__main__":
    main()
