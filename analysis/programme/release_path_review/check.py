"""Independent bounded smoke check of verify_bundle's path initialization."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import subprocess
from types import SimpleNamespace

BASELINE = "45ee5651dce6e668a17e39a956d8c1f673b41d09"


def check(repo):
    repo = Path(repo).resolve()
    source_path = repo / "analysis/verify_bundle.py"
    source = source_path.read_bytes()
    old = subprocess.check_output(["git", "show", BASELINE + ":analysis/verify_bundle.py"], cwd=repo)
    tree, original = ast.parse(source), ast.parse(old)
    scientific_ast = {}
    for name in ("require", "verify"):
        before = next(x for x in original.body if isinstance(x, ast.FunctionDef) and x.name == name)
        after = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == name)
        scientific_ast[name] = ast.dump(before, include_attributes=False) == ast.dump(after, include_attributes=False)
    before_main = next(x for x in original.body if isinstance(x, ast.If))
    after_main = next(x for x in tree.body if isinstance(x, ast.If) and ast.dump(x.test) == ast.dump(before_main.test))
    scientific_ast["main_entrypoint"] = ast.dump(before_main, include_attributes=False) == ast.dump(after_main, include_attributes=False)
    if not all(scientific_ast.values()):
        raise AssertionError("Scientific verification or main entrypoint changed")
    root_assignment = next(x for x in tree.body if isinstance(x, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "ROOT" for t in x.targets))
    guard = next(x for x in tree.body if isinstance(x, ast.If) and isinstance(x.test, ast.BoolOp))
    prefix_code = compile(ast.Module(body=[guard], type_ignores=[]), "<source-prefix-guard>", "exec")
    cases = []
    for name, initial, platform_name, expected in (
        ("drive_letter", "C:\\repo", "nt", "\\\\?\\C:\\repo"),
        ("already_extended_drive", "\\\\?\\C:\\repo", "nt", "\\\\?\\C:\\repo"),
        ("non_windows_guard", "C:\\repo", "posix", "C:\\repo"),
        ("unc", "\\\\server\\share\\repo", "nt", "\\\\?\\UNC\\server\\share\\repo"),
        ("already_extended_unc", "\\\\?\\UNC\\server\\share\\repo", "nt", "\\\\?\\UNC\\server\\share\\repo"),
    ):
        state = {"ROOT": PureWindowsPath(initial), "Path": PureWindowsPath, "os": SimpleNamespace(name=platform_name)}
        exec(prefix_code, state)
        actual = str(state["ROOT"])
        cases.append({"case": name, "initial": initial, "os_name": platform_name,
                      "actual": actual, "expected": expected, "passed": actual == expected})
    state = {"Path": Path, "__file__": str(source_path), "os": os}
    exec(compile(ast.Module(body=[root_assignment, guard], type_ignores=[]), "<source-root-initialization>", "exec"), state)
    effective = state["ROOT"]
    # This is the exact invocation used by the original audit. It must continue
    # to enumerate the same complete tracked inventory under the new root.
    plain_tree = subprocess.check_output(["git", "ls-tree", "-rz", "HEAD"], cwd=repo)
    extended_tree = subprocess.check_output(["git", "ls-tree", "-rz", "HEAD"], cwd=effective)
    if plain_tree != extended_tree:
        raise AssertionError("Tracked inventory changed under prefixed cwd")
    entries = []
    ordinary_missing, extended_missing = [], []
    for item in filter(None, plain_tree.split(b"\0")):
        header, encoded_name = item.split(b"\t", 1)
        name = encoded_name.decode("utf-8")
        expected = header.split()[2].decode("ascii")
        ordinary = repo / name
        corrected = effective / name
        if not ordinary.is_file():
            ordinary_missing.append(name)
        if not corrected.is_file():
            extended_missing.append(name)
        entries.append((len(str(ordinary)), name, expected))
    length, longest, expected_blob = max(entries)
    raw = (effective / longest).read_bytes()
    actual_blob = hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()
    if extended_missing or actual_blob != expected_blob:
        raise AssertionError("Corrected root cannot access all files or longest blob differs")
    if source_path.read_bytes() != source:
        raise AssertionError("Source changed during smoke check")
    return {
        "status": "passed_for_current_drive_letter_checkout" if all(x["passed"] for x in cases) else "current_checkout_passes_with_unc_portability_defect",
        "baseline_revision": BASELINE,
        "current_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo).decode().strip(),
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "baseline_source_sha256": hashlib.sha256(old).hexdigest(),
        "unchanged_ast": scientific_ast,
        "prefix_cases": cases,
        "same_complete_git_inventory": True,
        "tracked_count": len(entries),
        "ordinary_missing_count": len(ordinary_missing),
        "extended_missing": extended_missing,
        "longest_tracked_file": longest,
        "longest_ordinary_path_length": length,
        "longest_ordinary_is_file": (repo / longest).is_file(),
        "longest_extended_is_file": (effective / longest).is_file(),
        "longest_expected_git_blob": expected_blob,
        "longest_actual_git_blob": actual_blob,
        "scientific_verification_executed": False,
        "source_unchanged_during_review": True,
        "limits": "Smoke test and source review only. Does not replace the complete required integrity run after commit. UNC behavior tested lexically without network access.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    output = Path(args.out)
    if output.exists():
        raise FileExistsError("Fresh receipt path required")
    result = check(args.repo)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
