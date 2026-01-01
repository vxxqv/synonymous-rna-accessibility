"""Synthetic-only tests of the separately authored evaluator.

No observed feature, response, prediction or validation artifact is opened.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("independent_persist_evaluator", HERE / "evaluator.py")
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)
PROTOCOL = None
OUT = None
COUNTS = {"dual_ridge_equivalence_fits": 0, "independent_dual_cv_prediction_arrays": 0,
          "influence_records_checked": 0, "actual_outcome_or_feature_artifacts_read": 0}


def dual_ridge(x, y, future, penalty):
    """Separate small n-by-n solve, manually calculated training moments."""
    x, y, future = np.asarray(x, float), np.asarray(y, float), np.asarray(future, float)
    kept = [j for j in range(x.shape[1]) if len(set(map(float, x[:, j]))) > 1]
    means = [math.fsum(map(float, x[:, j])) / len(x) for j in kept]
    scales = [math.sqrt(math.fsum((float(v) - mean) ** 2 for v in x[:, j]) / len(x)) for j, mean in zip(kept, means)]
    z = np.asarray([[(float(row[j]) - mean) / scale for j, mean, scale in zip(kept, means, scales)] for row in x])
    w = np.asarray([[(float(row[j]) - mean) / scale for j, mean, scale in zip(kept, means, scales)] for row in future])
    center = math.fsum(map(float, y)) / len(y)
    alpha = np.linalg.solve(z @ z.T + penalty * np.eye(len(y)), y - center)
    return center + w @ z.T @ alpha


def fixture():
    ids = PROTOCOL["all_ids"]
    rng = np.random.default_rng(20261025)
    x = rng.normal(size=(37, 83))
    x[:, 0] = 7
    x[:, 5] = x[:, 4]
    x[:, 30] = 2 * x[:, 8]
    # Exact training constant, arbitrarily different on held-out rows.
    x[:22, 1] = 3
    x[22:, 1] = np.arange(15) + 1000
    y = {key: float((j * j + 3 * j) % 19 - 8) for j, key in enumerate(PROTOCOL["development_ids"])}
    values = {key: x[j].tolist() for j, key in enumerate(ids)}
    cds = {key: "A" * 666 for key in ids}
    return x, values, y, cds


class SyntheticTests(unittest.TestCase):
    def test_01_analytic_ridge(self):
        model = e.fit_augmented_svd([[-1], [1]], [-2, 2], 2)
        self.assertAlmostEqual(model["intercept"], 0, places=12)
        self.assertAlmostEqual(model["coefficients"][0], 1, places=12)
        np.testing.assert_allclose(e.predict(model, [[-1], [1]]), [-1, 1], atol=1e-12)

    def test_02_dual_equivalence_rank_deficiency(self):
        rng = np.random.default_rng(1367)
        for n, p in ((8, 21), (16, 82), (22, 83), (4, 1)):
            x = rng.normal(size=(n, p))
            future = rng.normal(size=(5, p))
            if p > 3:
                x[:, 0] = 9
                x[:, 2] = x[:, 1]
            y = rng.normal(size=n)
            for penalty in e.PENALTIES:
                fitted = e.fit_augmented_svd(x, y, penalty)
                expected = dual_ridge(x, y, future, penalty)
                np.testing.assert_allclose(e.predict(fitted, future), expected, rtol=1e-10, atol=1e-10)
                COUNTS["dual_ridge_equivalence_fits"] += 1

    def test_03_exact_rank_ties_and_constants(self):
        self.assertEqual(e.average_ranks([1, 2, 2, 4]), [1, 2.5, 2.5, 4])
        self.assertAlmostEqual(e.spearman([1, 2, 2, 4], [4, 2, 2, 1]), -1)
        self.assertIsNone(e.spearman([3, 3, 3], [1, 2, 3]))
        self.assertIsNone(e.spearman([1, 2, 3], [0, 0, 0]))
        self.assertEqual(e.average_ranks([1, 1 + 1e-14]), [1, 2])

    def test_04_prescribed_lambda_ties(self):
        y = [0.0, 10.0]
        tolerance = 1e-11
        losses = {0.1: 1.0, 1.0: 1.0 + tolerance / 4,
                  10.0: 1.0 + tolerance / 2, 100.0: 1.0 + 2 * tolerance}
        selected, _ = e.select_penalty(losses, y)
        self.assertEqual(selected, 10)
        scaled, _ = e.select_penalty({k: v * 200 for k, v in losses.items()}, [v * 200 for v in y])
        self.assertEqual(scaled, selected)
        exact, tie = e.select_penalty({0.1: 0, 1.0: 0, 10.0: 1e-16, 100.0: 0}, [0, 0])
        self.assertEqual(exact, 100)
        self.assertEqual(tie["absolute_loss_tolerance"], 0)

    def test_05_training_only_preparation(self):
        x = [[1, 10], [1, 11], [1, 12]]
        model = e.fit_augmented_svd(x, [1, 3, 4], 1)
        before = json.dumps(model, sort_keys=True)
        e.predict(model, [[100000, -1000000], [-100000, 1000000]])
        self.assertEqual(json.dumps(model, sort_keys=True), before)
        self.assertEqual(model["active_indices"], [1])
        self.assertEqual(model["means"], [11])
        self.assertAlmostEqual(model["scales"][0], math.sqrt(2 / 3))

    def test_06_full_fixed_cv_against_independent_dual(self):
        x, values, ymap, cds = fixture()
        result = e.fit_six(values, ymap, PROTOCOL, cds)
        ids = PROTOCOL["development_ids"]
        y = np.asarray([ymap[key] for key in ids])
        for name, width in e.RIDGE_WIDTHS.items():
            for penalty in e.PENALTIES:
                expected = np.zeros(22)
                for group in PROTOCOL["development_groups"]:
                    test = [j for j, key in enumerate(ids) if key in group]
                    train = [j for j, key in enumerate(ids) if key not in group]
                    expected[test] = dual_ridge(x[train, :width], y[train], x[test, :width], penalty)
                np.testing.assert_allclose(result["tuning"][name]["oof_predictions"][str(penalty)], expected, rtol=1e-10, atol=1e-10)
                self.assertAlmostEqual(result["tuning"][name]["mae_by_penalty"][str(penalty)], float(np.mean(abs(expected - y))), places=11)
                COUNTS["independent_dual_cv_prediction_arrays"] += 1
            selected = result["tuning"][name]["selected_penalty"]
            expected = dual_ridge(x[:22, :width], y, x[22:, :width], selected)
            np.testing.assert_allclose(result["predictions"][name], expected, rtol=1e-10, atol=1e-10)
        for key in PROTOCOL["holdout_ids"]:
            self.assertEqual(result["nearest_neighbors"][key], {"distance": 0, "training_ids": ids})
        np.testing.assert_allclose(result["predictions"]["nearest_CDS"], np.mean(y), rtol=0, atol=1e-14)
        # A deliberately extreme test-only replacement cannot alter CV or fit.
        changed = {key: list(row) for key, row in values.items()}
        for key in PROTOCOL["holdout_ids"]:
            changed[key] = [float(v * 1000 + 9999) for v in changed[key]]
        other = e.fit_six(changed, ymap, PROTOCOL, cds)
        self.assertEqual(result["fitted_ridges"], other["fitted_ridges"])
        self.assertEqual(result["tuning"], other["tuning"])

    def test_07_degenerate_response(self):
        _, values, _, cds = fixture()
        ymap = {key: 0.1 for key in PROTOCOL["development_ids"]}
        result = e.fit_six(values, ymap, PROTOCOL, cds)
        for name in e.RIDGE_WIDTHS:
            self.assertEqual(result["tuning"][name]["selected_penalty"], 100)
            self.assertEqual(set(result["tuning"][name]["mae_by_penalty"].values()), {0})
        for prediction in result["predictions"].values():
            np.testing.assert_allclose(prediction, 0.1, rtol=0, atol=1e-16)

    def test_08_zero_baseline_and_undefined_ranks(self):
        ids = PROTOCOL["holdout_ids"]
        record = {"holdout_ids": ids, "predictions": {name: [0.0] * 15 for name in e.MODELS}}
        result = e.evaluate(record, {key: 0.0 for key in ids})
        self.assertFalse(result["at_least_ten_percent_better_than_every_baseline"])
        self.assertEqual(result["zero_error_baselines"], list(e.BASELINES))
        self.assertTrue(all(value is None for value in result["fractional_gains"].values()))
        self.assertTrue(all(value["spearman_rank_agreement"] is None for value in result["models"].values()))
        self.assertTrue(all(value["undefined_omissions"] == 15 for value in result["influence_gain_ranges"].values()))

    def test_09_all_fifteen_influence_records(self):
        ids = PROTOCOL["holdout_ids"]
        targets = {key: float(j - 7) for j, key in enumerate(ids)}
        record = {"holdout_ids": ids, "predictions": {name: [float(j % (k + 2)) for j in range(15)] for k, name in enumerate(e.MODELS)}}
        result = e.evaluate(record, targets)
        y = np.asarray(list(targets.values()))
        self.assertEqual(len(result["leave_one_holdout_row_out"]), 15)
        for omitted, influence in enumerate(result["leave_one_holdout_row_out"]):
            mae = {name: float(np.mean(np.delete(abs(np.asarray(record["predictions"][name]) - y), omitted))) for name in e.MODELS}
            for name in e.MODELS:
                self.assertAlmostEqual(influence["mae"][name], mae[name], places=12)
            COUNTS["influence_records_checked"] += 1
        # One omission can produce an undefined gain even when the full gain is defined.
        record["predictions"] = {name: [0.0] * 15 for name in e.MODELS}
        for name in e.BASELINES:
            record["predictions"][name][0] = 1.0
        result = e.evaluate(record, {key: 0.0 for key in ids})
        self.assertTrue(result["at_least_ten_percent_better_than_every_baseline"])
        self.assertFalse(result["leave_one_holdout_row_out"][0]["at_least_ten_percent_better_than_every_baseline"])
        self.assertTrue(all(row["undefined_omissions"] == 1 for row in result["influence_gain_ranges"].values()))

    def test_10_invalid_groups_and_coverage_rejected(self):
        x, values, y, cds = fixture()
        with self.assertRaises(ValueError):
            e.tune(x[:22, :21], list(y.values()), PROTOCOL["development_ids"], [[PROTOCOL["development_ids"][0]]] * 22)
        incomplete = dict(y)
        incomplete.pop(next(iter(incomplete)))
        with self.assertRaises(ValueError):
            e.fit_six(values, incomplete, PROTOCOL, cds)
        with self.assertRaises(ValueError):
            e.fit_augmented_svd([[1], [2]], [1, float("nan")], 1)

    def test_11_authorization_precedes_observed_decoding(self):
        with tempfile.TemporaryDirectory(dir=OUT) as directory:
            root = Path(directory)
            auth = root / "authorization.json"
            auth.write_text(json.dumps({"status": "not_authorized", "operation": "independent_PERSIST_prediction_v2_compare"}))
            original_read = e.read
            opened = []
            def observed(path):
                opened.append(Path(path).name)
                return original_read(path)
            e.read = observed
            try:
                with self.assertRaises(ValueError):
                    e.verify_authorization(auth, root / "not_present_protocol.json")
            finally:
                e.read = original_read
            self.assertEqual(opened, ["authorization.json"])

    def test_12_dimensionless_metrics_do_not_use_response_units(self):
        with self.assertRaises(ValueError):
            e.compare_numeric({"spearman_rank_agreement": 0.5}, {"spearman_rank_agreement": 0.4}, 1e20)
        with self.assertRaises(ValueError):
            e.compare_numeric({"fractional_gains": {"baseline": 0.1}}, {"fractional_gains": {"baseline": 0.2}}, 1e20)
        with self.assertRaises(ValueError):
            e.compare_numeric({"influence_gain_ranges": {"baseline": {"minimum": 0.1}}},
                              {"influence_gain_ranges": {"baseline": {"minimum": 0.5}}}, 1e20)


def main():
    global PROTOCOL, OUT
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Do not overwrite synthetic evidence")
    if e.sha(args.protocol) != e.PROTOCOL_SHA256:
        raise ValueError("Wrong allowed protocol")
    PROTOCOL = e.read(args.protocol)
    args.out.mkdir(parents=True)
    OUT = args.out
    with (OUT / "unittest.txt").open("w", encoding="utf-8") as stream:
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SyntheticTests))
    receipt = {"status": "passed" if result.wasSuccessful() else "failed", "tests_run": result.testsRun,
               "failures": len(result.failures), "errors": len(result.errors), "counts": COUNTS,
               "inputs_sha256": {"protocol.json": e.sha(args.protocol), "evaluator.py": e.sha(HERE / "evaluator.py"),
                                   "test_synthetic.py": e.sha(Path(__file__))},
               "software": {"Python": sys.version.split()[0], "NumPy": np.__version__},
               "randomness": {"synthetic_only": True, "generator": "numpy.random.default_rng", "seeds": [1367, 20261025]},
               "independence": "No author implementation or observed numeric artifact imported/read; dual solves and explicit fixtures check the SVD evaluator"}
    e.write(OUT / "receipt.json", receipt)
    print(json.dumps({"status": receipt["status"], "tests_run": result.testsRun, "counts": COUNTS}))
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
