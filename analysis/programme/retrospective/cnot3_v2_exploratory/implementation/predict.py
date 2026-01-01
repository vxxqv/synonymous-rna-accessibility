"""Post-exposure exploratory CNOT3 operational-cohort comparison.

No workbook reader is provided. Training receives only training targets.
The CLI requires explicit hashes and authorization; see README.md.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import re
from typing import Any

import numpy as np

VERSION = "cnot3-operational-exploratory-v2"
LAMBDAS = (0.1, 1.0, 10.0, 100.0, 1000.0)
TIE_TOL = 1e-12
SENSE_CODONS = tuple("".join(c) for c in itertools.product("ACGT", repeat=3)
                     if "".join(c) not in ("TAA", "TAG", "TGA"))
AA = tuple("ACDEFGHIKLMNPQRSTVWY")
PRIMARY = ("no_change", "baseline_linear", "baseline_core", "amino_acid_core",
           "codon_core", "association")
MODELS = PRIMARY + ("codon_noSE", "association_noSE", "association_nonlinear")
SENSITIVITY_PAIRS = (("association_noSE", "codon_noSE"),
                     ("association_nonlinear", "codon_core"))
RAW_NUMERIC = ("HL_WT", "LOG2_FC_IP", "SE_IP", "CDS_length_nt", "CDS_GC",
               "UTR5_length_nt", "UTR5_GC", "UTR3_length_nt", "UTR3_GC")
FEATURE_KEYS = frozenset(("GeneID", "group_ID") + RAW_NUMERIC +
                         tuple("codon_" + c for c in SENSE_CODONS) +
                         tuple("aa_" + a for a in AA))
MINIMUMS = dict(training_genes=2000, holdout_genes=300, training_groups=100,
                holdout_groups=30, assessment_genes=100, inner_training_groups=50)


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def object_sha(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode("utf-8")).hexdigest()


def read_json(path: str | Path) -> Any:
    def unique_object(pairs: list[tuple]) -> dict:
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    def invalid_constant(value: str) -> None:
        raise ValueError(f"Nonfinite JSON constant: {value}")
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=unique_object, parse_constant=invalid_constant)


def write_json(path: str | Path, value: Any) -> None:
    with Path(path).open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def group_fold(group: str, partition: str) -> int:
    if partition not in ("outer", "inner"):
        raise ValueError("Unknown partition")
    encoded = ("cnot3-focus-v1|" + partition + "|" + group).encode("utf-8")
    return int.from_bytes(hashlib.sha256(encoded).digest(), "big") % 5


def finite(value: Any, label: str) -> float:
    if isinstance(value, (bool, str)) or value is None:
        raise ValueError(f"{label} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{label} must be a finite number") from error
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite {label}")
    return number


def validate_features(rows: list[dict], name: str) -> list[dict]:
    """Whitelist means a target hidden under an unexpected key is also rejected."""
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{name} must be a nonempty list")
    seen = set()
    output = []
    for original in rows:
        if not isinstance(original, dict) or set(original) != FEATURE_KEYS:
            extra = set(original) - FEATURE_KEYS if isinstance(original, dict) else set()
            missing = FEATURE_KEYS - set(original) if isinstance(original, dict) else FEATURE_KEYS
            raise ValueError(f"{name} feature schema differs: extra={sorted(extra)} missing={sorted(missing)}")
        row = dict(original)
        gene, group = row["GeneID"], row["group_ID"]
        if not isinstance(gene, str) or not gene or not isinstance(group, str) or not group:
            raise ValueError("GeneID and group_ID must be nonempty strings")
        if gene in seen:
            raise ValueError(f"Duplicate GeneID: {gene}")
        seen.add(gene)
        for key in FEATURE_KEYS - {"GeneID", "group_ID", "UTR5_GC", "UTR3_GC"}:
            row[key] = finite(row[key], key)
        for key in ("HL_WT", "SE_IP", "CDS_length_nt"):
            if row[key] <= 0:
                raise ValueError(f"Nonpositive {key}")
        for prefix in ("UTR5", "UTR3"):
            length = row[prefix + "_length_nt"]
            if length < 0:
                raise ValueError("Negative UTR length")
            # A zero-length UTR has no defined GC. Do not use any supplied value.
            if length == 0:
                if row[prefix + "_GC"] is not None:
                    gc = finite(row[prefix + "_GC"], prefix + "_GC")
                    if not 0 <= gc <= 1:
                        raise ValueError("GC outside [0,1]")
                row[prefix + "_GC"] = None
            else:
                row[prefix + "_GC"] = finite(row[prefix + "_GC"], prefix + "_GC")
                if not 0 <= row[prefix + "_GC"] <= 1:
                    raise ValueError("GC outside [0,1]")
        if not 0 <= row["CDS_GC"] <= 1:
            raise ValueError("CDS_GC outside [0,1]")
        for prefix, suffixes in (("codon_", SENSE_CODONS), ("aa_", AA)):
            values = [row[prefix + x] for x in suffixes]
            if any(not 0 <= x <= 1 for x in values):
                raise ValueError("Composition fraction outside [0,1]")
            # The reference schema excludes terminal stop from sense-codon fractions.
            # Do not introduce a sum constraint here; eligibility is source audited.
        output.append(row)
    return sorted(output, key=lambda x: x["GeneID"])


def ordered_targets(targets: dict[str, float], rows: list[dict]) -> np.ndarray:
    if not isinstance(targets, dict) or set(targets) != {r["GeneID"] for r in rows}:
        raise ValueError("Target GeneIDs must match the corresponding partition exactly")
    y = np.asarray([finite(targets[r["GeneID"]], "log(HL_KO)") for r in rows])
    # These are logs of positive finite half-lives, not arbitrary overflowing reals.
    with np.errstate(over="ignore", under="ignore"):
        raw = np.exp(y)
    if np.any(~np.isfinite(raw)) or np.any(raw <= 0):
        raise ValueError("Targets do not represent positive finite half-lives")
    return y


def check_design(train: list[dict], holdout: list[dict], receipt: dict) -> dict:
    train_ids = {r["GeneID"] for r in train}
    test_ids = {r["GeneID"] for r in holdout}
    train_groups = {r["group_ID"] for r in train}
    test_groups = {r["group_ID"] for r in holdout}
    if train_ids & test_ids or train_groups & test_groups:
        raise ValueError("Gene or group leakage across outer split")
    if any(group_fold(g, "outer") == 0 for g in train_groups):
        raise ValueError("Training group violates frozen hash partition")
    if any(group_fold(g, "outer") != 0 for g in test_groups):
        raise ValueError("Holdout group violates frozen hash partition")
    folds = np.asarray([group_fold(r["group_ID"], "inner") for r in train])
    counts = dict(training_genes=len(train), holdout_genes=len(holdout),
                  training_groups=len(train_groups), holdout_groups=len(test_groups))
    fold_counts = []
    for fold in range(5):
        assessment = int(np.sum(folds == fold))
        training_groups = len({r["group_ID"] for r, f in zip(train, folds) if f != fold})
        if assessment == 0 or training_groups == 0:
            raise ValueError("An inner fold has no assessment genes or training groups")
        fold_counts.append(dict(fold=fold, assessment_genes=assessment,
                                inner_training_groups=training_groups))
    synthetic = receipt.get("status") == "synthetic_only"
    if not synthetic:
        if receipt.get("status") != "root_authorized_operational_prediction":
            raise ValueError("No root authorization for observed training")
        if receipt.get("analysis_status") != "post_exposure_exploratory" or receipt.get("historical_published_cohort_reproduced") is not False or receipt.get("original_frozen_gate_status") != "held_for_source_diagnosis" or receipt.get("reproduction_gate_passed") is True:
            raise ValueError("Exploratory status and failed original reproduction must remain explicit")
        for field in ("operational_source_gate_sha256", "original_frozen_gate_sha256", "independent_source_audit_sha256"):
            if not re.fullmatch(r"[0-9a-f]{64}", str(receipt.get(field, ""))):
                raise ValueError("Missing hashed operational source evidence: " + field)
        for field in ("operational_source_gate_passed", "reference_eligibility_passed", "grouping_complete"):
            if receipt.get(field) is not True:
                raise ValueError(f"Required prefit gate is not passed: {field}")
        denominator = receipt.get("operational_slam_overlap_count")
        if not isinstance(denominator, int) or denominator < len(train) + len(holdout):
            raise ValueError("Invalid reference-coverage denominator")
        coverage = (len(train) + len(holdout)) / denominator
        if coverage < 0.8:
            raise ValueError("Reference coverage below 80 percent")
        if receipt.get("eligible_count") != len(train) + len(holdout):
            raise ValueError("Eligibility receipt count mismatch")
        for key, value in counts.items():
            if value < MINIMUMS[key]:
                raise ValueError(f"Insufficient {key}")
        for row in fold_counts:
            for key in ("assessment_genes", "inner_training_groups"):
                if row[key] < MINIMUMS[key]:
                    raise ValueError(f"Insufficient fold {key}")
    return dict(status=receipt["status"], counts=counts, fold_counts=fold_counts,
                source_receipt=receipt, source_receipt_sha256=object_sha(receipt))


def learn_transform(rows: list[dict], model: str) -> dict:
    x = np.log([r["HL_WT"] for r in rows])
    settings = dict(model=model, baseline_mean=float(np.mean(x)),
                    baseline_knots=np.quantile(x, [.1, .5, .9], method="linear").tolist())
    for prefix in ("UTR5", "UTR3"):
        values = [r[prefix + "_GC"] for r in rows if r[prefix + "_length_nt"] > 0]
        if not values:
            raise ValueError(f"No nonempty {prefix} in this training subset")
        settings[prefix + "_GC_median"] = float(np.median(values))
    settings["association_knots"] = (np.quantile([r["LOG2_FC_IP"] for r in rows],
                                              [.1, .5, .9], method="linear").tolist()
                                      if model == "association_nonlinear" else [])
    if model in ("association", "association_noSE", "association_nonlinear"):
        base, _ = base_matrix(rows, settings)
        mean, scale = base.mean(axis=0), base.std(axis=0, ddof=0)
        constant = np.all(base == base[0], axis=0)
        scale[constant] = 1.0
        if np.any(scale <= 0) or not np.all(np.isfinite(scale)):
            raise ValueError("Invalid projection-base scale")
        z = (base - mean) / scale
        z[:, constant] = 0.0
        augmented = np.column_stack((np.ones(len(rows)), z))
        raw, names = association_basis(rows, settings)
        projections = []
        for column, name in zip(raw.T, names):
            raw_mean = float(column.mean())
            raw_constant = bool(np.all(column == column[0]))
            raw_scale = 1.0 if raw_constant else float(column.std(ddof=0))
            standardized = np.zeros(len(rows)) if raw_constant else (column - raw_mean) / raw_scale
            coefficients, _, rank, singular = np.linalg.lstsq(augmented, standardized, rcond=1e-12)
            residual = standardized - augmented @ coefficients
            norm = float(np.linalg.norm(residual))
            zero = bool(norm <= 1e-10 * math.sqrt(len(rows)))
            residual_mean = float(residual.mean())
            residual_scale = float(residual.std(ddof=0))
            if not zero and (residual_scale <= 0 or not math.isfinite(residual_scale)):
                raise ValueError("Unusable association residual scale")
            projections.append(dict(name=name, raw_mean=raw_mean, raw_scale=raw_scale,
                                    raw_constant=raw_constant, coefficients=coefficients.tolist(),
                                    rank=int(rank), singular_values=singular.tolist(),
                                    residual_norm=norm, residual_mean=residual_mean,
                                    residual_scale=residual_scale, zero=zero))
        settings["association_projection"] = dict(base_mean=mean.tolist(), base_scale=scale.tolist(),
                                                  base_constant=constant.tolist(), columns=projections,
                                                  rcond=1e-12, cutoff_per_sqrt_n=1e-10)
    return settings


def base_matrix(rows: list[dict], settings: dict) -> tuple[np.ndarray, list[str]]:
    model = settings["model"]
    columns, names = [], []
    def add(name: str, values: Any) -> None:
        names.append(name)
        columns.append(np.asarray(values, dtype=float))
    baseline = np.log([r["HL_WT"] for r in rows])
    add("log_HL_WT", baseline)
    centered = baseline - settings["baseline_mean"]
    add("centered_log_HL_WT_squared", centered ** 2)
    add("centered_log_HL_WT_cubed", centered ** 3)
    for index, knot in enumerate(settings["baseline_knots"]):
        add(f"baseline_hinge_{index}", np.maximum(baseline - knot, 0))
    add("log_CDS_length", np.log([r["CDS_length_nt"] for r in rows]))
    add("CDS_GC", [r["CDS_GC"] for r in rows])
    for prefix in ("UTR5", "UTR3"):
        add("log1p_" + prefix + "_length", np.log1p([r[prefix + "_length_nt"] for r in rows]))
        add(prefix + "_GC", [r[prefix + "_GC"] if r[prefix + "_length_nt"] > 0
                              else settings[prefix + "_GC_median"] for r in rows])
    for prefix in ("UTR5", "UTR3"):
        add(prefix + "_empty", [float(r[prefix + "_length_nt"] == 0) for r in rows])
    if model not in ("codon_noSE", "association_noSE"):
        add("log_SE_IP", np.log([r["SE_IP"] for r in rows]))
    if model == "amino_acid_core":
        for aa in AA:
            add("aa_" + aa, [r["aa_" + aa] for r in rows])
    if model in ("codon_core", "association", "codon_noSE", "association_noSE", "association_nonlinear"):
        for codon in SENSE_CODONS:
            add("codon_" + codon, [r["codon_" + codon] for r in rows])
    matrix = np.column_stack(columns)
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Nonfinite transformed feature")
    return matrix, names


def association_basis(rows: list[dict], settings: dict) -> tuple[np.ndarray, list[str]]:
    association = np.asarray([r["LOG2_FC_IP"] for r in rows])
    columns, names = [association], ["LOG2_FC_IP"]
    for index, knot in enumerate(settings["association_knots"]):
        columns.append(np.maximum(association - knot, 0))
        names.append(f"association_hinge_{index}")
    return np.column_stack(columns), names


def design_matrix(rows: list[dict], settings: dict) -> tuple[np.ndarray, list[str]]:
    base, names = base_matrix(rows, settings)
    if "association_projection" not in settings:
        return base, names
    projection = settings["association_projection"]
    z = (base - np.asarray(projection["base_mean"])) / np.asarray(projection["base_scale"])
    z[:, np.asarray(projection["base_constant"], dtype=bool)] = 0.0
    augmented = np.column_stack((np.ones(len(rows)), z))
    raw, _ = association_basis(rows, settings)
    residuals = []
    for column, fitted in zip(raw.T, projection["columns"]):
        if fitted["zero"]:
            # This applies even if an assessment row breaks the training duplicate.
            residuals.append(np.zeros(len(rows)))
        else:
            standardized = (column - fitted["raw_mean"]) / fitted["raw_scale"]
            residual = standardized - augmented @ np.asarray(fitted["coefficients"])
            residuals.append((residual - fitted["residual_mean"]) / fitted["residual_scale"])
        names.append("residualized_" + fitted["name"])
    matrix = np.column_stack((base, *residuals))
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Nonfinite residualized feature")
    return matrix, names


def ridge_prepare(x: np.ndarray, y: np.ndarray) -> dict:
    mean = x.mean(axis=0)
    scale = x.std(axis=0, ddof=0)
    constant = np.all(x == x[0], axis=0)
    # Exact constant columns are zero even when floating reduction creates tiny SD.
    scale[constant] = 1.0
    if np.any(scale <= 0) or not np.all(np.isfinite(scale)):
        raise ValueError("Invalid training scale")
    z = (x - mean) / scale
    z[:, constant] = 0.0
    u, singular, vt = np.linalg.svd(z, full_matrices=False)
    return dict(mean=mean, scale=scale, constant=constant, y_mean=float(np.mean(y)),
                singular=singular, vt=vt, uy=u.T @ (y - np.mean(y)))


def ridge_parameters(prepared: dict, lam: float) -> dict:
    singular = prepared["singular"]
    beta = prepared["vt"].T @ (singular / (singular ** 2 + lam) * prepared["uy"])
    return dict(mean=prepared["mean"].tolist(), scale=prepared["scale"].tolist(),
                constant=prepared["constant"].tolist(), y_mean=prepared["y_mean"],
                beta=beta.tolist(), penalty=float(lam))


def ridge_predict(x: np.ndarray, fit: dict) -> np.ndarray:
    z = (x - np.asarray(fit["mean"])) / np.asarray(fit["scale"])
    z[:, np.asarray(fit["constant"], dtype=bool)] = 0.0
    predicted = fit["y_mean"] + z @ np.asarray(fit["beta"])
    if not np.all(np.isfinite(predicted)):
        raise ValueError("Nonfinite prediction")
    return predicted


def choose_lambda(scores: dict[float, float]) -> float:
    if set(scores) != set(LAMBDAS) or not all(math.isfinite(x) for x in scores.values()):
        raise ValueError("Invalid tuning scores")
    best = min(scores.values())
    return max(lam for lam in LAMBDAS if scores[lam] <= best + TIE_TOL)


def fit_ridge_model(train: list[dict], y: np.ndarray, holdout: list[dict], model: str) -> dict:
    folds = np.asarray([group_fold(r["group_ID"], "inner") for r in train])
    oof = {lam: np.full(len(train), np.nan) for lam in LAMBDAS}
    details = []
    for fold in range(5):
        fit_mask, assess_mask = folds != fold, folds == fold
        fit_rows = [r for r, keep in zip(train, fit_mask) if keep]
        assess_rows = [r for r, keep in zip(train, assess_mask) if keep]
        transform = learn_transform(fit_rows, model)
        x, names = design_matrix(fit_rows, transform)
        assessment_x, _ = design_matrix(assess_rows, transform)
        prepared = ridge_prepare(x, y[fit_mask])
        fits = {}
        for lam in LAMBDAS:
            fit = ridge_parameters(prepared, lam)
            oof[lam][assess_mask] = ridge_predict(assessment_x, fit)
            fits[str(lam)] = fit
        details.append(dict(fold=fold, transform=transform, columns=names,
                            training_ids=[r["GeneID"] for r in fit_rows],
                            assessment_ids=[r["GeneID"] for r in assess_rows], fits=fits))
    if not all(np.all(np.isfinite(v)) for v in oof.values()):
        raise ValueError("Incomplete OOF prediction")
    scores = {lam: float(np.mean(np.abs(pred - y))) for lam, pred in oof.items()}
    chosen = choose_lambda(scores)
    transform = learn_transform(train, model)
    x, names = design_matrix(train, transform)
    holdout_x, _ = design_matrix(holdout, transform)
    fitted = ridge_parameters(ridge_prepare(x, y), chosen)
    return dict(selected_lambda=chosen, pooled_oof_mae={str(k): v for k, v in scores.items()},
                oof_predictions={str(k): v.tolist() for k, v in oof.items()},
                fold_fits=details, transform=transform, columns=names, fitted=fitted,
                holdout_prediction=ridge_predict(holdout_x, fitted).tolist())


def predict(training_features: list[dict], training_targets: dict[str, float],
            holdout_features: list[dict], design_receipt: dict) -> dict:
    """No argument can carry held-out targets. Unknown predictor fields fail closed."""
    train = validate_features(training_features, "training")
    holdout = validate_features(holdout_features, "holdout")
    y = ordered_targets(training_targets, train)
    design = check_design(train, holdout, design_receipt)
    x = np.log([r["HL_WT"] for r in train])
    holdout_x = np.log([r["HL_WT"] for r in holdout])
    linear_mean = float(x.mean())
    linear_constant = bool(np.all(x == x[0]))
    linear_scale = 1.0 if linear_constant else float(x.std(ddof=0))
    z = np.zeros(len(x)) if linear_constant else (x - linear_mean) / linear_scale
    assessment_z = np.zeros(len(holdout_x)) if linear_constant else (holdout_x - linear_mean) / linear_scale
    coefficients, _, rank, singular = np.linalg.lstsq(np.column_stack((np.ones(len(z)), z)), y, rcond=1e-12)
    models = {
        "no_change": dict(holdout_prediction=holdout_x.tolist()),
        "baseline_linear": dict(training_mean=linear_mean, training_scale=linear_scale,
                                constant=linear_constant, coefficients=coefficients.tolist(),
                                rank=int(rank), singular_values=singular.tolist(), rcond=1e-12,
                                holdout_prediction=(coefficients[0] + coefficients[1] * assessment_z).tolist()),
    }
    for model in MODELS[2:]:
        models[model] = fit_ridge_model(train, y, holdout, model)
    # Store training responses for the separately specified whole-cohort description.
    # Nothing here contains a held-out response.
    context = lambda rows: [dict(GeneID=r["GeneID"], group_ID=r["group_ID"],
                                HL_WT=r["HL_WT"], LOG2_FC_IP=r["LOG2_FC_IP"]) for r in rows]
    return dict(version=VERSION, status=design["status"], model_order=list(MODELS),
                runtime=dict(python=platform.python_version(), numpy=np.__version__),
                design=design, primary_baselines=list(PRIMARY[:-1]),
                training_ids=[r["GeneID"] for r in train],
                holdout_ids=[r["GeneID"] for r in holdout],
                training_fold_ids=[group_fold(r["group_ID"], "inner") for r in train],
                training_context=context(train), holdout_context=context(holdout),
                training_logKO=y.tolist(),
                baseline_quartiles=np.quantile(x, [.25, .5, .75], method="linear").tolist(),
                h_ref=float(np.median([r["HL_WT"] for r in train])), models=models,
                input_sha256=dict(training_features=object_sha(training_features),
                                  training_targets=object_sha(training_targets),
                                  holdout_features=object_sha(holdout_features)))


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + j - 1) / 2 + 1
        i = j
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 2 or len(y) != len(x):
        return None
    a, b = average_ranks(np.asarray(x)), average_ranks(np.asarray(y))
    a, b = a - a.mean(), b - b.mean()
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.clip(np.dot(a, b) / denom, -1, 1)) if denom else None


def metrics(pred: np.ndarray, y: np.ndarray, groups: list[str]) -> dict:
    signed = pred - y
    absolute = np.abs(signed)
    return dict(MAE=float(np.mean(absolute)), RMSE=float(np.sqrt(np.mean(signed ** 2))),
                median_absolute_error=float(np.median(absolute)),
                p90_absolute_error=float(np.quantile(absolute, .9, method="linear")),
                mean_signed_error=float(np.mean(signed)), spearman=spearman(pred, y),
                equal_component_MAE=float(np.mean([np.mean(absolute[np.asarray(groups) == g])
                                                   for g in sorted(set(groups))])))


def relative_gain(primary: float, baseline: float) -> float | None:
    return (baseline - primary) / baseline if baseline > 0 else None


def comparisons(errors: dict[str, float]) -> dict:
    primary = errors["association"]
    result = {}
    for model in PRIMARY[:-1]:
        base = errors[model]
        result[model] = dict(absolute_MAE_reduction=base - primary,
                             relative_MAE_reduction=relative_gain(primary, base),
                             meets_5_percent=bool(base > 0 and primary <= 0.95 * base))
    return result


def kinetic_summary(context: list[dict], y: np.ndarray, h_ref: float) -> dict:
    wt = np.asarray([r["HL_WT"] for r in context])
    ip = np.asarray([r["LOG2_FC_IP"] for r in context])
    logwt = np.log(wt)
    log2change = (y - logwt) / math.log(2)
    # Algebraic direct definition. Fail transparently on unrepresentable values.
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        rate = math.log(2) * (h_ref / wt - h_ref * np.exp(-y))
    if not np.all(np.isfinite(rate)):
        raise ValueError("Nonfinite dimensionless rate contrast")
    arrays = dict(baseline_HL=wt, log2_lifetime_change=log2change,
                  dimensionless_rate_contrast=rate, LOG2_FC_IP=ip)
    summaries = {}
    for name, values in arrays.items():
        summaries[name] = dict(median=float(np.median(values)),
                               p10=float(np.quantile(values, .1, method="linear")),
                               p90=float(np.quantile(values, .9, method="linear")),
                               negative=int(np.sum(values < 0)), zero=int(np.sum(values == 0)),
                               positive=int(np.sum(values > 0)))
    correlations = {a + "__" + b: spearman(arrays[a], arrays[b])
                    for a, b in itertools.combinations(arrays, 2)}
    return dict(n=len(y), h_ref=h_ref, summaries=summaries, spearman=correlations,
                rows=[dict(GeneID=row["GeneID"], **{key: float(values[i]) for key, values in arrays.items()})
                      for i, row in enumerate(context)])


def evaluate(predictions: dict, holdout_targets: dict[str, float]) -> dict:
    """Pure arithmetic evaluator. CLI enforces the external prediction-commit receipt."""
    if predictions.get("version") != VERSION or tuple(predictions.get("model_order", [])) != MODELS:
        raise ValueError("Prediction model/version mismatch")
    context = predictions["holdout_context"]
    if predictions["holdout_ids"] != [r["GeneID"] for r in context]:
        raise ValueError("Holdout ordering mismatch")
    y = ordered_targets(holdout_targets, context)
    groups = [r["group_ID"] for r in context]
    fields = {}
    predicted = {}
    for model in MODELS:
        pred = np.asarray(predictions["models"][model]["holdout_prediction"], dtype=float)
        if pred.shape != y.shape or not np.all(np.isfinite(pred)):
            raise ValueError("Malformed frozen predictions")
        predicted[model] = pred
        fields[model] = metrics(pred, y, groups)
    comparison = comparisons({model: fields[model]["MAE"] for model in PRIMARY})
    success = all(row["meets_5_percent"] for row in comparison.values())
    quartiles = np.asarray(predictions["baseline_quartiles"])
    bins = np.searchsorted(quartiles, np.log([r["HL_WT"] for r in context]), side="left")
    strata = []
    for interval in range(4):
        mask = bins == interval
        strata.append(dict(interval=interval, n=int(np.sum(mask)),
                           lower=None if interval == 0 else float(quartiles[interval - 1]),
                           upper=None if interval == 3 else float(quartiles[interval]),
                           convention="lower open, upper closed",
                           MAE={m: float(np.mean(np.abs(p[mask] - y[mask]))) if np.any(mask) else None
                                for m, p in predicted.items()}))
    deletions = []
    for group in sorted(set(groups)):
        mask = np.asarray(groups) != group
        if not np.any(mask):
            raise ValueError("Cannot delete sole holdout component")
        errors = {model: float(np.mean(np.abs(predicted[model][mask] - y[mask]))) for model in PRIMARY}
        compared = comparisons(errors)
        deletions.append(dict(deleted_group=group, remaining_genes=int(np.sum(mask)),
                              MAE=errors, comparisons=compared,
                              primary_success=all(r["meets_5_percent"] for r in compared.values())))
    gains = [d["comparisons"]["codon_core"]["relative_MAE_reduction"] for d in deletions]
    absolute = [d["comparisons"]["codon_core"]["absolute_MAE_reduction"] for d in deletions]
    sensitivities = {}
    for model, baseline in SENSITIVITY_PAIRS:
        base, primary = fields[baseline]["MAE"], fields[model]["MAE"]
        sensitivities[model] = dict(baseline=baseline, absolute_MAE_reduction=base - primary,
                                    relative_MAE_reduction=relative_gain(primary, base))
    train_context = predictions["training_context"]
    train_y = np.asarray(predictions["training_logKO"])
    all_context, all_y = train_context + context, np.concatenate((train_y, y))
    rows = [dict(GeneID=row["GeneID"], group_ID=row["group_ID"], log_HL_KO=float(y[i]),
                 baseline_interval=int(bins[i]),
                 models={m: dict(predicted=float(p[i]), signed_error=float(p[i] - y[i]),
                                  absolute_error=float(abs(p[i] - y[i]))) for m, p in predicted.items()})
            for i, row in enumerate(context)]
    return dict(version=VERSION, status=predictions["status"], primary_success=success,
                interpretation="Model-held-out processed-measurement prediction within the named comparators; not causal information or independent replication",
                metrics=fields, primary_comparisons=comparison, sensitivities=sensitivities,
                baseline_strata=strata, rows=rows, component_deletions=deletions,
                influence=dict(codon_absolute_gain_range=[min(absolute), max(absolute)],
                               codon_relative_gain_range=([min(gains), max(gains)] if all(x is not None for x in gains) else None),
                               undefined_relative_deletions=sum(x is None for x in gains),
                               positive_codon_gain_every_deletion=all(x > 0 for x in absolute),
                               positive_gain_every_baseline_every_deletion=all(
                                   row["absolute_MAE_reduction"] > 0 for d in deletions for row in d["comparisons"].values()),
                               five_percent_every_deletion=all(d["primary_success"] for d in deletions)),
                kinetic=dict(entire_eligible=kinetic_summary(all_context, all_y, predictions["h_ref"]),
                             holdout=kinetic_summary(context, y, predictions["h_ref"])))


def verify_input_freeze(paths: dict[str, str | Path], freeze_path: str | Path, statuses: tuple[str, ...]) -> dict:
    freeze = read_json(freeze_path)
    if freeze.get("status") not in statuses:
        raise ValueError("Freeze does not authorize this phase")
    if freeze.get("implementation_sha256") != sha256(__file__):
        raise ValueError("Implementation hash mismatch")
    for label, path in paths.items():
        if freeze.get("inputs", {}).get(label) != sha256(path):
            raise ValueError(f"Frozen input hash mismatch: {label}")
    return freeze


def run_predict(args: argparse.Namespace) -> None:
    out = Path(args.out)
    if out.exists():
        raise FileExistsError("Output directory already exists")
    paths = dict(training_features=args.training_features, training_targets=args.training_targets,
                 holdout_features=args.holdout_features, design_receipt=args.design_receipt,
                 protocol=args.protocol)
    freeze = verify_input_freeze(paths, args.freeze, ("synthetic_only", "root_authorized_operational_prediction"))
    receipt = read_json(args.design_receipt)
    if receipt.get("status") != freeze["status"]:
        raise ValueError("Design and freeze status mismatch")
    out.mkdir(parents=True, exist_ok=False)
    try:
        result = predict(read_json(args.training_features), read_json(args.training_targets),
                         read_json(args.holdout_features), receipt)
        result["file_provenance"] = dict(implementation_sha256=sha256(__file__),
                                          freeze_sha256=sha256(args.freeze), inputs=freeze["inputs"])
        write_json(out / "predictions.json", result)
        write_json(out / "receipt.json", dict(status="predictions_saved_unscored", origin=freeze["status"],
                                              predictions_sha256=sha256(out / "predictions.json"),
                                              implementation_sha256=sha256(__file__)))
    except Exception as error:
        write_json(out / "failure.json", dict(error_type=type(error).__name__, message=str(error)))
        raise


def run_evaluate(args: argparse.Namespace) -> None:
    out = Path(args.out)
    if out.exists():
        raise FileExistsError("Output directory already exists")
    # Hash and commit assertions precede loading any held-out target file.
    actual = sha256(args.predictions)
    if actual != args.expected_predictions_sha256:
        raise ValueError("Prediction hash mismatch before target access")
    commit = read_json(args.commit_receipt)
    if commit.get("predictions_sha256") != actual:
        raise ValueError("Commit receipt prediction hash mismatch")
    predictions = read_json(args.predictions)
    synthetic = predictions.get("status") == "synthetic_only"
    if synthetic:
        if commit.get("status") != "synthetic_only":
            raise ValueError("Synthetic evaluation requires synthetic receipt")
    elif (commit.get("status") != "root_verified_operational_prediction_commit" or
          not re.fullmatch(r"[0-9a-f]{40}", str(commit.get("commit", ""))) or
          commit.get("verified") is not True):
        raise ValueError("Root must verify exact prediction commit before evaluation")
    if predictions.get("file_provenance", {}).get("implementation_sha256") != sha256(__file__):
        raise ValueError("Evaluation implementation differs from prediction freeze")
    if commit.get("holdout_targets_sha256") != sha256(args.holdout_targets):
        raise ValueError("Held-out target file hash mismatch before decoding")
    out.mkdir(parents=True, exist_ok=False)
    try:
        result = evaluate(predictions, read_json(args.holdout_targets))
        result["provenance"] = dict(predictions_sha256=actual,
                                     holdout_targets_sha256=sha256(args.holdout_targets),
                                     commit_receipt_sha256=sha256(args.commit_receipt),
                                     implementation_sha256=sha256(__file__))
        write_json(out / "evaluation.json", result)
        write_json(out / "receipt.json", dict(status="evaluation_complete", origin=predictions["status"],
                                              evaluation_sha256=sha256(out / "evaluation.json")))
    except Exception as error:
        write_json(out / "failure.json", dict(error_type=type(error).__name__, message=str(error)))
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    train = commands.add_parser("predict", help="Fit using training targets only and save locked predictions")
    for option in ("training-features", "training-targets", "holdout-features", "design-receipt", "protocol", "freeze", "out"):
        train.add_argument("--" + option, required=True)
    score = commands.add_parser("evaluate", help="Score after the root-verified prediction commit")
    for option in ("predictions", "expected-predictions-sha256", "commit-receipt", "holdout-targets", "out"):
        score.add_argument("--" + option, required=True)
    args = parser.parse_args()
    (run_predict if args.phase == "predict" else run_evaluate)(args)


if __name__ == "__main__":
    main()
