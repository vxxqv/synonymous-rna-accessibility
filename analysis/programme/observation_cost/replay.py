"""Portable replay of the frozen 384-case observation-cost study.

Requires preinstalled NetworkX 3.7. No Git, installation or network calls occur.
Mathematical functions and the independent reference remain byte-preserved.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from fractions import Fraction as Q
import gzip
import hashlib
import importlib.util
import itertools
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
BASE_REL = Path("work/cost_audit_cycle2")
ACQUISITION_SHA = "6005d30d709dc8100a1de64822d8e80cf8a8dbc8d9c3a18de8d21f69622be6ee"
FREEZE_SHA = "3686acf9a6f5c1210cff47630de0959f536d1c552d57577f7157c939bae8a9d3"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def relative(name):
    path = PurePosixPath(name.replace("\\", "/"))
    require(bool(path.parts) and not path.is_absolute(), f"Not a relative path: {name!r}")
    require(all(p not in (".", "..") and ":" not in p for p in path.parts), f"Unsafe path: {name!r}")
    return Path(*path.parts)


def hashes(root, paths):
    return {p: sha(root / relative(p)) for p in paths}


def load_module(name, source):
    spec = importlib.util.spec_from_file_location(name, source)
    require(spec is not None and spec.loader is not None, "Could not import preserved source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def first_difference(a, b, at="$"):
    if type(a) is not type(b):
        return f"{at}: different types"
    if isinstance(a, dict):
        if a.keys() != b.keys():
            return f"{at}: different key sets"
        for key in a:
            difference = first_difference(a[key], b[key], f"{at}.{key}")
            if difference:
                return difference
    elif isinstance(a, list):
        if len(a) != len(b):
            return f"{at}: different lengths"
        for i, (left, right) in enumerate(zip(a, b)):
            difference = first_difference(left, right, f"{at}[{i}]")
            if difference:
                return difference
    elif a != b:
        return f"{at}: {a!r} != {b!r}"
    return None


def compare(expected, fresh, ignored=()):
    require(set(expected) == set(fresh), "Different receipt top-level keys")
    require(set(ignored) <= set(expected), "Unknown ignored receipt fields")
    left = {k: v for k, v in expected.items() if k not in ignored}
    right = {k: v for k, v in fresh.items() if k not in ignored}
    difference = first_difference(left, right)
    require(difference is None, f"Scientific mismatch: {difference}")
    return {"exact_agreement": True, "ignored_top_level_fields": list(ignored),
            "compared_top_level_fields": sorted(left)}


def replay_author(base):
    """Original case order and orchestration; all mathematics calls held code.

    The original run() wraps this body in Windows-path and Git provenance code.
    This adapter preserves its loop, checks, runtime cap and serializations.
    """
    module = load_module("portable_cost_author", base / "cost_audit.py")
    output = base / "results"
    require(not output.exists(), "Fresh author results must not already exist")
    output.mkdir()
    start = perf_counter()
    counts, rows, certificates = Counter(), [], []
    small = module.small_checks()
    module.dump(output / "c4_bruteforce.json", small)
    adjacency = tuple(tuple(j for j in range(8) if (i ^ j).bit_count() == 1) for i in range(8))
    model = module.GraphModel(adjacency, (0, 7), range(1, 7), ((1, 6), (2, 5), (3, 4)))
    require(len(model.states) == 729 and len(model.worlds) == 64, "Wrong frozen graph state counts")
    counts["partial_states_checked"] = len(model.states)
    counts["completion_certificate_checks"] = model.completion_checks
    protocol = read(base / "protocol.json")
    for costs, p, batch in itertools.product(itertools.product((1, 4), repeat=6),
                                            (Q(1, 4), Q(1, 2), Q(3, 4)), (False, True)):
        if perf_counter() - start > protocol["runtime_cap_seconds"]:
            raise TimeoutError("Original frozen author runtime cap reached")
        values, policy, checks = module.solve(model, costs, p, batch)
        counts["bellman_action_inequalities"] += checks
        counts["bellman_states"] += len(values)
        case_id = len(rows)
        row = {"case_id": case_id, "singleton_costs": costs, "p": str(p),
               "mode": "singleton_plus_pairs" if batch else "singleton",
               "optimal_expected_cost": str(values[model.root]), "policies": {}}
        for name in module.POLICIES:
            item = module.evaluate(model, costs, p, batch, name, policy)
            item["regret_to_optimum"] = str(Q(item["expected_cost"]) - values[model.root])
            require(Q(item["regret_to_optimum"]) >= 0, "Negative regret")
            if name == "exact_bayes_optimal":
                require(Q(item["expected_cost"]) == values[model.root], "Optimal replay mismatch")
            row["policies"][name] = item
            counts["full_world_policy_checks"] += len(model.worlds)
        rows.append(row)
        certificates.append({"case_id": case_id, "states": [
            [module.key(state), str(values[state]), model.action_names[policy[state]] if state in policy else None]
            for state in model.states]})
    require(len(rows) == 384 and counts["full_world_policy_checks"] == 384 * 64 * 4, "Incomplete frozen grid")
    summary = module.summarize(rows)
    summary["checks"] = dict(counts)
    module.dump(output / "cases.json", {"world_order": model.worlds,
                "world_result_columns": ["paid_cost", "ordered_action_names", "newly_revealed_bits"], "cases": rows})
    module.dump(output / "summary.json", summary)
    raw = json.dumps(certificates, separators=(",", ":")).encode()
    (output / "bellman_certificates.json.gz").write_bytes(gzip.compress(raw, mtime=0))
    return {"runtime_seconds": perf_counter() - start, "checks": dict(counts)}


def reference_worker(root):
    """Execute unchanged reference with a narrowly scoped path-only adapter.

    Only json.loads of the exact frozen freeze text normalizes its relative path
    keys in memory. Other JSON, including every scientific value, is untouched.
    """
    base = root / BASE_REL
    freeze_path = base / "freeze.json"
    require(sha(freeze_path) == FREEZE_SHA, "Wrong reference freeze")
    freeze_text = freeze_path.read_text(encoding="utf-8")
    original_loads = json.loads
    calls = []

    def portable_loads(text, *args, **kwargs):
        obj = original_loads(text, *args, **kwargs)
        if text == freeze_text:
            obj["input_sha256"] = {relative(p).as_posix(): h for p, h in obj["input_sha256"].items()}
            calls.append("only frozen metadata relative-path keys normalized")
        return obj

    json.loads = portable_loads
    try:
        load_module("portable_cost_independent", base / "reviewer/independent_reference.py")
    finally:
        json.loads = original_loads
    require(len(calls) == 1, "Reference metadata adapter did not run exactly once")
    write(base / "reviewer/path_adapter_receipt.json", {"status": "passed", "adaptations": calls,
          "source_bytes_changed": False, "scientific_JSON_changed": False})


def reaggregate_reference(base, preserved):
    """Repeat the reviewed exact aggregation from independent case expectations."""
    reference = read(base / "reviewer/independent_results.json")
    rows = read(base / "results/cases.json")["cases"]
    expected = read(base / "results/summary.json")
    values = reference["comparison_results"]
    require(len(rows) == len(values) == 384, "Wrong independent summary grid")
    require(all(a["case_id"] == b["case_id"] for a, b in zip(rows, values)), "Independent case order mismatch")
    policies = {}
    for name in expected["policies"]:
        costs = [Q(r["policy_expected_costs"][name]) for r in values]
        optimum = [Q(r["optimal_expected_cost"]) for r in values]
        regret = [a - b for a, b in zip(costs, optimum)]
        policies[name] = {"optimal_ties": regret.count(0), "strict_losses": sum(x > 0 for x in regret),
            "mean_regret_over_fixed_grid": str(sum(regret, Q(0)) / 384), "maximum_regret": str(max(regret)),
            "maximum_ratio": str(max(a / b for a, b in zip(costs, optimum)))}
    grouped = {}
    for row, independent in zip(rows, values):
        grouped.setdefault((tuple(row["singleton_costs"]), row["p"]), {})[row["mode"]] = Q(independent["optimal_expected_cost"])
    gain = [p["singleton"] - p["singleton_plus_pairs"] for p in grouped.values()]
    batch = {"pairs": len(gain), "strict_improvements": sum(x > 0 for x in gain), "ties": gain.count(0),
             "mean_gain_over_fixed_grid": str(sum(gain, Q(0)) / len(gain)), "maximum_gain": str(max(gain))}
    fresh = {"cases": 384, "policies": policies, "optional_batch_value": batch}
    for left in ("unchanged_path_cut_singletons", "cheapest_cost_per_unseen_bit"):
        right = "unchanged_lazy_forward_singletons"
        difference = [Q(r["policy_expected_costs"][left]) - Q(r["policy_expected_costs"][right]) for r in values]
        fresh[left + "_minus_" + right] = {"wins": sum(x < 0 for x in difference), "ties": difference.count(0),
            "losses": sum(x > 0 for x in difference), "minimum_difference": str(min(difference)),
            "maximum_difference": str(max(difference)), "mean_difference_over_fixed_grid": str(sum(difference, Q(0)) / 384)}
    compare({k: v for k, v in expected.items() if k not in ("checks", "scope")}, fresh)
    verification = {"status": "pass", "policies": policies, "optional_batch_value": batch,
                    "all_saved_summary_fields_recomputed": True}
    comparison = compare(read(preserved / BASE_REL / "reviewer/summary_verification.json"), verification)
    write(base / "reviewer/summary_verification.json", verification)
    return comparison


def run(out, manifest):
    preserved, active = out / "preserved", out / "reconstruction"
    preserved.mkdir()
    active.mkdir()
    with zipfile.ZipFile(HERE / "frozen_sources.zip") as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)) and set(names) == set(manifest["payload_sha256"]), "Wrong ZIP inventory")
        for info in archive.infolist():
            rel = relative(info.filename)
            require(not info.is_dir() and ((info.external_attr >> 16) & 0o170000) != 0o120000, "Invalid ZIP member")
            dest = preserved / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(archive.read(info))
    before = hashes(preserved, manifest["payload_sha256"])
    require(before == manifest["payload_sha256"], "Frozen payload differs")
    require(before["outputs/synonymous-rna-accessibility/src/acquisition.py"] == ACQUISITION_SHA, "Wrong selector source")
    frozen = read(preserved / BASE_REL / "freeze.json")
    for path, expected in frozen["input_sha256"].items():
        require(before[relative(path).as_posix()] == expected, f"Freeze mismatch: {path}")
    for name in manifest["active_input_paths"]:
        rel = relative(name)
        dest = active / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(preserved / rel, dest)
    active_before = hashes(active, manifest["active_input_paths"])
    base = active / BASE_REL
    with (out / "author_stdout.txt").open("w", encoding="utf-8") as stdout, \
            (out / "author_stderr.txt").open("w", encoding="utf-8") as stderr, \
            redirect_stdout(stdout), redirect_stderr(stderr):
        author_receipt = replay_author(base)
    comparisons = {}
    for name in ("cases.json", "summary.json", "c4_bruteforce.json"):
        comparisons[name] = compare(read(preserved / BASE_REL / "results" / name), read(base / "results" / name))
    original_cert = gzip.decompress((preserved / BASE_REL / "results/bellman_certificates.json.gz").read_bytes())
    fresh_cert = gzip.decompress((base / "results/bellman_certificates.json.gz").read_bytes())
    require(original_cert == fresh_cert, "Exact Bellman certificate serialization differs")
    comparisons["bellman_certificates.json.gz"] = {"exact_decompressed_bytes_agreement": True,
        "raw_gzip_bytes_agree": sha(preserved / BASE_REL / "results/bellman_certificates.json.gz") == sha(base / "results/bellman_certificates.json.gz"),
        "ignored_serialization_metadata": "Gzip wrapper compression/header bytes may vary by host; decompressed bytes must be identical"}
    reviewer = base / "reviewer"
    reviewer.mkdir(exist_ok=True)
    require(not (reviewer / "independent_results.json").exists(), "Historical reference receipt in active output")
    start = perf_counter()
    completed = subprocess.run([sys.executable, "-I", "-B", str(Path(__file__).resolve()), "--reference-root", str(active)],
                               cwd=active, capture_output=True, text=True, encoding="utf-8", timeout=240, check=False)
    reference_seconds = perf_counter() - start
    (out / "reference_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (out / "reference_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    require(completed.returncode == 0, f"Independent reference exited {completed.returncode}")
    fresh_reference = read(reviewer / "independent_results.json")
    comparisons["independent_results.json"] = compare(
        read(preserved / BASE_REL / "reviewer/independent_results.json"), fresh_reference,
        ("inputs_sha256", "command", "python", "runtime_seconds", "started_at"))
    # Input-hash provenance is not compared across different absolute paths and
    # host text serialization. It is independently checked against actual bytes.
    for path, expected in fresh_reference["inputs_sha256"].items():
        resolved = Path(path).resolve()
        require(resolved.is_relative_to(active), "Reference followed a path outside reconstructed inputs")
        require(sha(resolved) == expected, "Reference input hash is stale")
    require(fresh_reference["input_hashes_unchanged"] is True, "Reference input changed")
    comparisons["independent_summary_verification.json"] = reaggregate_reference(base, preserved)
    after = hashes(preserved, manifest["payload_sha256"])
    active_after = hashes(active, manifest["active_input_paths"])
    require(before == after and active_before == active_after, "Immutable input changed")
    require(not list(out.rglob("__pycache__")) and not list(out.rglob(".git")), "Unexpected cache or Git metadata")
    output_paths = [str((BASE_REL / "results" / name).as_posix()) for name in
                    ("cases.json", "summary.json", "c4_bruteforce.json", "bellman_certificates.json.gz")]
    output_paths += [(BASE_REL / "reviewer" / name).as_posix() for name in
                     ("independent_results.json", "summary_verification.json", "path_adapter_receipt.json")]
    return {"status": "passed", "scientific_comparisons": comparisons,
        "frozen_inputs_unchanged": True, "payload_sha256_before": before, "payload_sha256_after": after,
        "active_input_sha256_before": active_before, "active_input_sha256_after": active_after,
        "historical_freeze_input_count": len(frozen["input_sha256"]),
        "historical_git_verified": False, "git_required": False, "network_used": False,
        "runtime_seconds": {"author": author_receipt["runtime_seconds"], "independent_reference": reference_seconds},
        "creator_finite_checks": author_receipt["checks"], "independent_finite_checks": fresh_reference["counts"],
        "generated_result_sha256": hashes(active, output_paths),
        "scope": "384 frozen full-support cases only; no mathematical method, correlated-prior experiment or new branch added"}


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O or -OO")
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--out", type=Path, help="New output directory, which must not already exist")
    group.add_argument("--reference-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.reference_root is not None:
        reference_worker(args.reference_root.resolve())
        return 0
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError("Refusing to overwrite an existing output directory")
    manifest = read(HERE / "integration_manifest.json")
    for name, expected in manifest["package_sha256"].items():
        require(sha(HERE / relative(name)) == expected, f"Package mismatch: {name}")
    import networkx as nx
    require(nx.__version__ == "3.7", "Preinstalled NetworkX 3.7 is required; replay never installs dependencies")
    out.mkdir(parents=True)
    receipt = {"schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
               "environment": {"Python": platform.python_version(), "NetworkX": nx.__version__,
                               "hardware": platform.platform() + "; " + platform.machine()},
               "package_sha256": manifest["package_sha256"], "assertions_enabled": __debug__,
               "historical_git_note": "Original manifests retained as history; no historical Git event reauthenticated"}
    try:
        receipt.update(run(out, manifest))
    except Exception as error:
        receipt.update({"status": "failed", "error": f"{type(error).__name__}: {error}"})
        (out / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write(out / "receipt.json", receipt)
        print(json.dumps({"status": "failed", "error": receipt["error"]}))
        return 1
    write(out / "receipt.json", receipt)
    print(json.dumps({"status": "passed", "cases": 384, "exact_scientific_agreement": True,
                      "historical_git_verified": False, "receipt": "receipt.json"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
