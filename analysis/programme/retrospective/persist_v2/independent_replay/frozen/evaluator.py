"""Outcome-blind independent evaluator for the frozen PERSIST prediction protocol.

No author implementation is imported. Ridge uses an augmented SVD least-squares
problem with an unpenalized intercept. Observed inputs require an explicit hash
authorization, including this evaluator's own source hash.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import traceback

import numpy as np

sys.dont_write_bytecode = True
PROTOCOL_SHA256 = "a3dc686a23acaa12af2f2f44154b8baac819cb6df3a0a6d81328dd3ebdc5fd24"
AUTHOR_INTERFACE_SHA256 = "962f953ed17574db96f5ba0cde4827486830244941ea745bb34fb697a5514aaa"
MODELS = ("development_median", "development_mean", "composition_ridge",
          "codon_composition_ridge", "nearest_CDS", "codon_composition_MUP_ridge")
RIDGE_WIDTHS = {MODELS[2]: 21, MODELS[3]: 82, MODELS[5]: 83}
PENALTIES = (0.1, 1.0, 10.0, 100.0)
PRIMARY = MODELS[-1]
BASELINES = MODELS[:-1]
RELATIVE_COMPARISON_TOLERANCE = 1e-9
ABSOLUTE_RESPONSE_SCALE_FACTOR = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def array(values, ndim, label):
    result = np.asarray(values, dtype=float)
    require(result.ndim == ndim and np.all(np.isfinite(result)), f"Invalid finite {label}")
    return result


def fit_augmented_svd(x, y, penalty):
    x, y = array(x, 2, "training design"), array(y, 1, "training responses")
    require(len(x) == len(y) and len(y) > 0, "Training shape mismatch")
    require(math.isfinite(penalty) and penalty > 0, "Positive penalty required")
    # Exact equality determines column removal, not an estimated variance cutoff.
    indices = np.flatnonzero(np.any(x != x[0], axis=0))
    selected = x[:, indices]
    means = selected.mean(axis=0)
    scales = np.sqrt(np.mean((selected - means) ** 2, axis=0))
    require(np.all(scales > 0) and np.all(np.isfinite(scales)), "Invalid nonconstant-column scale")
    z = (selected - means) / scales
    p = len(indices)
    if np.all(y == y[0]):
        # This is the exact mathematical minimizer, including nonbinary decimal
        # constants whose repeated floating-point mean may round differently.
        intercept, coefficients = float(y[0]), np.zeros(p)
    elif p == 0:
        intercept, coefficients = math.fsum(map(float, y)) / len(y), np.zeros(0)
    else:
        design = np.zeros((len(y) + p, p + 1))
        design[:len(y), 0] = 1
        design[:len(y), 1:] = z
        design[len(y):, 1:] = math.sqrt(penalty) * np.eye(p)
        response = np.concatenate((y, np.zeros(p)))
        solution = np.linalg.lstsq(design, response, rcond=None)[0]
        intercept, coefficients = float(solution[0]), solution[1:]
    require(math.isfinite(intercept) and np.all(np.isfinite(coefficients)), "Nonfinite fitted model")
    return {"penalty": float(penalty), "active_indices": indices.tolist(), "means": means.tolist(),
            "scales": scales.tolist(), "coefficients": coefficients.tolist(), "intercept": intercept,
            "training_count": len(y), "feature_count": x.shape[1], "standardization_ddof": 0,
            "method": "augmented least squares, SVD, unpenalized intercept"}


def predict(model, x):
    x = array(x, 2, "prediction design")
    require(x.shape[1] == model["feature_count"], "Wrong prediction width")
    z = (x[:, model["active_indices"]] - np.asarray(model["means"])) / np.asarray(model["scales"])
    beta = model["coefficients"]
    result = np.asarray([math.fsum([model["intercept"], *(float(a) * b for a, b in zip(row, beta))]) for row in z])
    require(np.all(np.isfinite(result)), "Nonfinite prediction")
    return result


def select_penalty(losses, y):
    y = array(y, 1, "development responses")
    require(set(losses) == set(PENALTIES), "Incomplete penalty grid")
    require(all(math.isfinite(v) and v >= 0 for v in losses.values()), "Invalid CV loss")
    span = float(max(y) - min(y))
    tolerance = 0.0 if span == 0 else 1e-12 * span
    best = min(losses.values())
    eligible = [penalty for penalty in PENALTIES if losses[penalty] <= best + tolerance]
    return max(eligible), {"response_range": span, "absolute_loss_tolerance": tolerance,
                           "minimum_loss": best, "eligible_penalties": eligible}


def tune(x, y, ids, groups):
    x, y = array(x, 2, "CV design"), array(y, 1, "CV responses")
    require(len(ids) == len(y) == len(x) and len(set(ids)) == len(ids), "CV row mismatch")
    flattened = [key for group in groups for key in group]
    require(len(flattened) == len(ids) and set(flattened) == set(ids), "Groups must partition development IDs")
    require(all(group for group in groups) and len(groups) > 1, "Degenerate grouping")
    where = {key: j for j, key in enumerate(ids)}
    losses, oof, folds = {}, {}, {}
    for penalty in PENALTIES:
        predictions = np.full(len(ids), np.nan)
        preparations = []
        for group in groups:
            test = np.asarray([where[key] for key in group], dtype=int)
            training = np.asarray([j for j, key in enumerate(ids) if key not in set(group)], dtype=int)
            require(len(training) > 0, "Empty training fold")
            model = fit_augmented_svd(x[training], y[training], penalty)
            predictions[test] = predict(model, x[test])
            preparations.append({"omitted_ids": list(group), "training_ids": [ids[j] for j in training],
                                 "model": model})
        require(np.all(np.isfinite(predictions)), "Incomplete out-of-fold predictions")
        losses[penalty] = math.fsum(abs(float(a) - float(b)) for a, b in zip(predictions, y)) / len(y)
        oof[str(penalty)] = predictions.tolist()
        folds[str(penalty)] = preparations
    selected, tie = select_penalty(losses, y)
    return fit_augmented_svd(x, y, selected), {"selected_penalty": selected,
        "mae_by_penalty": {str(k): v for k, v in losses.items()}, "tie": tie,
        "oof_predictions": oof, "folds": folds}


def fit_six(values, development_y, protocol, cds):
    train, holdout = protocol["development_ids"], protocol["holdout_ids"]
    require(set(development_y) == set(train), "Require all and only development responses")
    require(set(values) == set(train + holdout) and set(cds) == set(values), "Complete fixed cohort required")
    require(tuple(protocol["models"]) == MODELS and tuple(protocol["ridge_lambdas"]) == PENALTIES, "Changed model contract")
    x_train = array([values[key] for key in train], 2, "development features")
    x_test = array([values[key] for key in holdout], 2, "holdout features")
    require(x_train.shape[1] == x_test.shape[1] == 83, "Require all 83 fixed columns")
    y = array([development_y[key] for key in train], 1, "development responses")
    predictions = {MODELS[0]: [float(np.median(y))] * len(holdout),
                   MODELS[1]: [math.fsum(map(float, y)) / len(y)] * len(holdout)}
    fitted, tuning = {}, {}
    for name, width in RIDGE_WIDTHS.items():
        fitted[name], tuning[name] = tune(x_train[:, :width], y, train, protocol["development_groups"])
        predictions[name] = predict(fitted[name], x_test[:, :width]).tolist()
    require(all(isinstance(cds[key], str) and cds[key] for key in cds), "Invalid CDS metadata")
    normalized = {key: cds[key].upper().replace("T", "U") for key in cds}
    require(len({len(sequence) for sequence in normalized.values()}) == 1, "Different CDS lengths")
    require(all(set(sequence) <= set("ACGU") for sequence in normalized.values()), "Unexpected CDS alphabet")
    neighbors, nearest = {}, []
    for key in holdout:
        distances = [sum(a != b for a, b in zip(normalized[key], normalized[t])) for t in train]
        minimum = min(distances)
        ties = [j for j, distance in enumerate(distances) if distance == minimum]
        nearest.append(math.fsum(float(y[j]) for j in ties) / len(ties))
        neighbors[key] = {"distance": minimum, "training_ids": [train[j] for j in ties]}
    predictions["nearest_CDS"] = nearest
    return {"development_ids": train, "holdout_ids": holdout, "models": list(MODELS),
            "predictions": {name: predictions[name] for name in MODELS}, "fitted_ridges": fitted,
            "tuning": tuning, "nearest_neighbors": neighbors,
            "development_targets": {key: float(development_y[key]) for key in train}}


def average_ranks(values):
    values = array(values, 1, "rank values")
    order = sorted(range(len(values)), key=lambda j: float(values[j]))
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        rank = (start + 1 + stop) / 2
        for j in order[start:stop]:
            result[j] = rank
        start = stop
    return result


def spearman(a, b):
    require(len(a) == len(b) and len(a) > 0, "Rank shape mismatch")
    ra, rb = average_ranks(a), average_ranks(b)
    ca = math.fsum(ra) / len(ra)
    cb = math.fsum(rb) / len(rb)
    da, db = [r - ca for r in ra], [r - cb for r in rb]
    aa, bb = math.fsum(d * d for d in da), math.fsum(d * d for d in db)
    if aa == 0 or bb == 0:
        return None
    return math.fsum(a * b for a, b in zip(da, db)) / math.sqrt(aa * bb)


def gains(maes):
    require(set(maes) == set(MODELS) and all(math.isfinite(v) and v >= 0 for v in maes.values()), "Invalid MAEs")
    zero = [name for name in BASELINES if maes[name] == 0]
    fractional = {name: None if maes[name] == 0 else (maes[name] - maes[PRIMARY]) / maes[name] for name in BASELINES}
    return {"at_least_ten_percent_better_than_every_baseline": not zero and all(maes[PRIMARY] <= 0.9 * maes[name] for name in BASELINES),
            "fractional_gains": fractional, "zero_error_baselines": zero}


def evaluate(record, targets):
    ids = record["holdout_ids"]
    require(set(targets) == set(ids) and len(ids) >= 2, "Complete holdout targets required")
    y = array([targets[key] for key in ids], 1, "holdout responses")
    summaries, errors = {}, {}
    for name in MODELS:
        prediction = array(record["predictions"][name], 1, "fixed predictions")
        require(len(prediction) == len(y), "Wrong prediction count")
        delta = [float(a) - float(b) for a, b in zip(prediction, y)]
        errors[name] = [abs(d) for d in delta]
        summaries[name] = {"mae": math.fsum(errors[name]) / len(y),
            "rmse": math.sqrt(math.fsum(d * d for d in delta) / len(y)),
            "spearman_rank_agreement": spearman(prediction, y)}
    influence = []
    for omitted, key in enumerate(ids):
        mae = {name: math.fsum(v for j, v in enumerate(errors[name]) if j != omitted) / (len(y) - 1) for name in MODELS}
        influence.append({"omitted_id": key, "mae": mae, **gains(mae)})
    ranges = {}
    for name in BASELINES:
        defined = [row["fractional_gains"][name] for row in influence if row["fractional_gains"][name] is not None]
        ranges[name] = {"minimum": min(defined) if defined else None, "maximum": max(defined) if defined else None,
                        "undefined_omissions": len(influence) - len(defined)}
    return {"models": summaries, **gains({name: summaries[name]["mae"] for name in MODELS}),
        "rows": [{"SequenceID": key, "target": float(y[j]),
                  "predictions": {name: record["predictions"][name][j] for name in MODELS},
                  "absolute_errors": {name: errors[name][j] for name in MODELS}} for j, key in enumerate(ids)],
        "leave_one_holdout_row_out": influence, "influence_gain_ranges": ranges,
        "influence_interpretation": "Dependent sensitivity diagnostics, not confidence intervals",
        "estimand": "Deposited cellular degradation coefficient, original numeric units"}


def compare_numeric(expected, actual, absolute_tolerance, at="$", stats=None):
    if stats is None:
        stats = {"numeric_leaves": 0, "exact_leaves": 0, "maximum_absolute_numeric_difference": 0.0}
    if expected is None or isinstance(expected, (str, bool)):
        require(expected == actual, f"Exact mismatch at {at}")
        stats["exact_leaves"] += 1
    elif isinstance(expected, (float, int)):
        require(isinstance(actual, (float, int)) and not isinstance(actual, bool), f"Wrong numeric type at {at}")
        require(math.isfinite(expected) and math.isfinite(actual), f"Nonfinite comparison at {at}")
        difference = abs(float(expected) - float(actual))
        dimensionless = any(field in at for field in ("fractional_gains", "influence_gain_ranges", "spearman_rank_agreement"))
        local_absolute = 1e-12 if dimensionless else absolute_tolerance
        local_relative = 1e-10 if dimensionless else RELATIVE_COMPARISON_TOLERANCE
        require(difference <= local_absolute + local_relative * abs(float(expected)), f"Numeric mismatch at {at}")
        stats["numeric_leaves"] += 1
        stats["maximum_absolute_numeric_difference"] = max(stats["maximum_absolute_numeric_difference"], difference)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(expected) == len(actual), f"List mismatch at {at}")
        for j, (a, b) in enumerate(zip(expected, actual)):
            compare_numeric(a, b, absolute_tolerance, f"{at}[{j}]", stats)
    elif isinstance(expected, dict):
        require(isinstance(actual, dict) and expected.keys() == actual.keys(), f"Dictionary mismatch at {at}")
        for key in expected:
            compare_numeric(expected[key], actual[key], absolute_tolerance, f"{at}.{key}", stats)
    else:
        raise ValueError(f"Unsupported comparison at {at}")
    return stats


def verify_authorization(path, protocol_path):
    authorization = read(path)
    require(authorization.get("status") == "authorized_by_root" and authorization.get("operation") == "independent_PERSIST_prediction_v2_compare", "No observed-data authorization")
    require(authorization.get("independent_evaluator_sha256") == sha(__file__), "Evaluator not frozen in authorization")
    require(sha(protocol_path) == PROTOCOL_SHA256 == authorization.get("protocol_sha256"), "Wrong frozen protocol")
    required = {"features", "cds", "selected_targets", "author_predictions", "author_validation"}
    require(set(authorization.get("artifacts", {})) == required, "Explicit artifact authorization required")
    verified = {}
    for name, artifact in authorization["artifacts"].items():
        p = Path(artifact["path"]).resolve()
        require(p.is_file() and sha(p) == artifact["sha256"], f"Authorized artifact hash mismatch: {name}")
        verified[name] = p
    return authorization, verified


def run_comparison(protocol_path, authorization_path, out):
    authorization, paths = verify_authorization(authorization_path, protocol_path)
    before = {name: sha(path) for name, path in paths.items()}
    protocol = read(protocol_path)
    features = read(paths["features"])
    require(features["protocol_sha256"] == PROTOCOL_SHA256 and features["ids"] == protocol["all_ids"] and features["names"] == protocol["feature_names"], "Wrong feature identity or order")
    cds_record = read(paths["cds"])
    require(set(cds_record) == {"cds_by_id"}, "Require CDS-only metadata")
    targets = read(paths["selected_targets"])
    require(set(targets) == {"endpoint", "development_targets", "holdout_targets"} and targets["endpoint"] == protocol["endpoint"], "Require selected-target-only artifact")
    require(set(targets["development_targets"]) == set(protocol["development_ids"]) and set(targets["holdout_targets"]) == set(protocol["holdout_ids"]), "Target coverage mismatch")
    # Do not decode author numeric answers until independent reconstruction is done.
    prediction = fit_six(features["values"], targets["development_targets"], protocol, cds_record["cds_by_id"])
    evaluation = evaluate(prediction, targets["holdout_targets"])
    write(out / "independent_predictions.json", prediction)
    write(out / "independent_validation.json", evaluation)
    author_prediction = read(paths["author_predictions"])
    author_evaluation = read(paths["author_validation"])
    require(author_prediction["protocol_sha256"] == PROTOCOL_SHA256 and author_prediction["features_sha256"] == before["features"]
            and author_prediction["implementation_sha256"] == AUTHOR_INTERFACE_SHA256
            and author_prediction["holdout_targets_read"] is False, "Author prediction provenance or exposure mismatch")
    require(author_evaluation["protocol_sha256"] == PROTOCOL_SHA256, "Wrong validation protocol")
    require(set(author_evaluation) == set(evaluation) | {"prediction_freeze", "protocol_sha256"}, "Unexpected validation schema")
    require(prediction["nearest_neighbors"] == author_prediction["nearest_neighbors"], "Nearest-CDS tie/distance mismatch")
    scale = max(1.0, *(abs(float(v)) for v in targets["development_targets"].values()),
                *(abs(float(v)) for v in targets["holdout_targets"].values()))
    absolute_tolerance = ABSOLUTE_RESPONSE_SCALE_FACTOR * scale
    results = {}
    for field in ("development_ids", "holdout_ids", "models", "predictions", "nearest_neighbors", "development_targets"):
        results[field] = compare_numeric(prediction[field], author_prediction[field], absolute_tolerance, field)
    for name in RIDGE_WIDTHS:
        expected, observed = prediction["tuning"][name], author_prediction["tuning"][name]
        require(expected["selected_penalty"] == observed["selected_penalty"], f"Penalty-selection mismatch for {name}")
        for field in ("mae_by_penalty", "oof_predictions"):
            results[name + "." + field] = compare_numeric(expected[field], observed[field], absolute_tolerance, name + "." + field)
    results["validation_all_scientific_fields"] = compare_numeric(
        evaluation, {key: author_evaluation[key] for key in evaluation}, absolute_tolerance, "validation")
    after = {name: sha(path) for name, path in paths.items()}
    require(before == after, "Observed input changed during independent comparison")
    return {"status": "passed", "protocol_sha256": PROTOCOL_SHA256, "evaluator_sha256": sha(__file__),
            "author_interface_sha256": AUTHOR_INTERFACE_SHA256, "artifacts_sha256_before": before,
            "artifacts_sha256_after": after, "all_inputs_unchanged": True, "comparisons": results,
            "numeric_tolerance": {"absolute": absolute_tolerance, "relative": RELATIVE_COMPARISON_TOLERANCE,
                 "scale_formula": "max(1, absolute development/holdout responses)",
                 "dimensionless_metrics": {"absolute": 1e-12, "relative": 1e-10},
                 "exact_requirements": ["selected penalty", "Boolean decisions", "None/undefined status", "IDs and grouping/order"]},
            "independence_limit": "MUP and fixed feature values are specified frozen inputs; no refolding or independent structural-feature provenance claim",
            "scope": "Fixed-cohort numerical reproduction, not calibration, causal structure, biological independence or a new predictive endpoint"}


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    require(not args.out.exists(), "Refusing to overwrite output")
    # Authorization and all file hashes are verified before creating an output
    # directory or decoding any observed numeric artifact.
    verify_authorization(args.authorization, args.protocol)
    args.out.mkdir(parents=True)
    receipt = {"started_at": datetime.now(timezone.utc).isoformat(),
               "software": {"Python": platform.python_version(), "NumPy": np.__version__}}
    try:
        receipt.update(run_comparison(args.protocol, args.authorization, args.out))
    except Exception as error:
        receipt.update({"status": "failed", "error": f"{type(error).__name__}: {error}", "retuned": False})
        (args.out / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        write(args.out / "receipt.json", receipt)
        return 1
    write(args.out / "receipt.json", receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
