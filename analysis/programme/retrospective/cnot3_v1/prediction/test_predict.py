"""Synthetic tests only. This file never opens the source biological workbooks."""
import argparse
import copy
import inspect
import itertools
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import predict as p

CODE = dict(zip(("".join(c) for c in itertools.product("TCAG", repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def synthetic(seed=741, n_groups=180, members=2):
    rng = np.random.default_rng(seed)
    train, test, train_y, test_y = [], [], {}, {}
    for group_index in range(n_groups):
        group = f"SYN{group_index:05d}_0"
        for member in range(members):
            gene = f"SYN{group_index:05d}_{member}"
            x = rng.uniform(-1.5, 1.5)
            ip = rng.normal()
            codons = rng.dirichlet(np.full(61, 3.0))
            aas = {aa: sum(float(v) for codon, v in zip(p.SENSE_CODONS, codons) if CODE[codon] == aa)
                   for aa in p.AA}
            sense_length = int(rng.integers(100, 1200))
            gc = sum(float(v) * (codon.count("G") + codon.count("C"))
                     for codon, v in zip(p.SENSE_CODONS, codons)) * sense_length / (3 * (sense_length + 1))
            row = dict(GeneID=gene, group_ID=group, HL_WT=math.exp(x),
                       LOG2_FC_IP=ip, SE_IP=math.exp(rng.uniform(-1, 1)),
                       CDS_length_nt=float((sense_length + 1) * 3), CDS_GC=gc,
                       UTR5_length_nt=float(rng.integers(10, 800)), UTR5_GC=float(rng.uniform(.2, .8)),
                       UTR3_length_nt=float(rng.integers(10, 1200)), UTR3_GC=float(rng.uniform(.2, .8)))
            if (group_index + member) % 7 == 0:
                row.update(UTR5_length_nt=0.0, UTR5_GC=None)
            if (group_index + member) % 11 == 0:
                row.update(UTR3_length_nt=0.0, UTR3_GC=None)
            row.update({"codon_" + c: float(v) for c, v in zip(p.SENSE_CODONS, codons)})
            row.update({"aa_" + a: aas[a] for a in p.AA})
            # A deliberately informative independent association, plus baseline cubic.
            y = x + .15 * x ** 2 - .08 * x ** 3 + 1.7 * ip
            if p.group_fold(group, "outer") == 0:
                test.append(row)
                test_y[gene] = y
            else:
                train.append(row)
                train_y[gene] = y
    return train, train_y, test, test_y


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.train, cls.y, cls.holdout, cls.hy = synthetic()
        cls.result = p.predict(cls.train, cls.y, cls.holdout, {"status": "synthetic_only"})
        cls.evaluation = p.evaluate(cls.result, cls.hy)

    def test_complete_models_and_informative_signal(self):
        self.assertEqual(tuple(self.result["models"]), p.MODELS)
        self.assertEqual(len(self.evaluation["rows"]), len(self.holdout))
        self.assertTrue(self.evaluation["primary_success"])
        self.assertGreater(self.evaluation["primary_comparisons"]["codon_core"]["relative_MAE_reduction"], .5)
        self.assertEqual(set(self.evaluation["sensitivities"]), {"association_noSE", "association_nonlinear"})

    def test_no_heldout_targets_in_training_api_or_predictors(self):
        self.assertNotIn("holdout_targets", inspect.signature(p.predict).parameters)
        bad = copy.deepcopy(self.holdout)
        bad[0]["HL_CNOT3KO"] = 1.0
        with mock.patch.object(p, "fit_ridge_model", side_effect=AssertionError("fit must not begin")):
            with self.assertRaisesRegex(ValueError, "extra=.*HL_CNOT3KO"):
                p.predict(self.train, self.y, bad, {"status": "synthetic_only"})
        y = dict(self.y, **self.hy)
        with self.assertRaisesRegex(ValueError, "Target GeneIDs"):
            p.predict(self.train, y, self.holdout, {"status": "synthetic_only"})

    def test_fold_preprocessing_and_projection_use_only_training(self):
        rows = p.validate_features(self.train, "synthetic")
        for model in p.MODELS[2:]:
            for fold in self.result["models"][model]["fold_fits"]:
                fit = [r for r in rows if p.group_fold(r["group_ID"], "inner") != fold["fold"]]
                expected = p.learn_transform(fit, model)
                self.assertEqual(expected, fold["transform"])
                self.assertFalse(set(fold["training_ids"]) & set(fold["assessment_ids"]))
                fit_groups = {r["group_ID"] for r in rows if r["GeneID"] in fold["training_ids"]}
                assess_groups = {r["group_ID"] for r in rows if r["GeneID"] in fold["assessment_ids"]}
                self.assertFalse(fit_groups & assess_groups)

    def test_assessment_extremes_do_not_change_fit(self):
        train = p.validate_features(self.train, "train")
        test = p.validate_features(self.holdout, "holdout")
        y = p.ordered_targets(self.y, train)
        altered = copy.deepcopy(test)
        for row in altered:
            row["HL_WT"] = 1e6
            row["LOG2_FC_IP"] = 50.0
            row["UTR5_GC"] = .999 if row["UTR5_length_nt"] else None
        other = p.fit_ridge_model(train, y, altered, "association")
        original = self.result["models"]["association"]
        for key in ("selected_lambda", "pooled_oof_mae", "oof_predictions", "fold_fits", "transform", "fitted"):
            self.assertEqual(original[key], other[key])
        self.assertNotEqual(original["holdout_prediction"], other["holdout_prediction"])

    def test_utr_imputation_training_only_and_missing_design(self):
        rows = p.validate_features(self.train, "train")
        settings = p.learn_transform(rows, "codon_core")
        expected = np.median([r["UTR5_GC"] for r in rows if r["UTR5_length_nt"] > 0])
        self.assertEqual(settings["UTR5_GC_median"], expected)
        missing = copy.deepcopy(rows)
        for row in missing:
            row["UTR5_length_nt"], row["UTR5_GC"] = 0.0, None
        with self.assertRaisesRegex(ValueError, "No nonempty UTR5"):
            p.learn_transform(missing, "codon_core")
        invalid = copy.deepcopy(self.train)
        invalid[0]["CDS_GC"] = None
        with self.assertRaises(ValueError):
            p.validate_features(invalid, "bad")

    def test_duplicate_and_near_duplicate_cutoff(self):
        rows = p.validate_features(self.train, "train")
        for epsilon in (0.0, 1e-13):
            duplicate = copy.deepcopy(rows)
            for i, row in enumerate(duplicate):
                row["LOG2_FC_IP"] = row["CDS_GC"] + epsilon * ((i % 2) * 2 - 1)
            settings = p.learn_transform(duplicate, "association")
            projection = settings["association_projection"]["columns"][0]
            self.assertTrue(projection["zero"])
            x, _ = p.design_matrix(duplicate, settings)
            np.testing.assert_array_equal(x[:, -1], 0.0)
            changed = copy.deepcopy(duplicate[:5])
            for row in changed:
                row["LOG2_FC_IP"] = 1000.0
            assessed, _ = p.design_matrix(changed, settings)
            np.testing.assert_array_equal(assessed[:, -1], 0.0)
        nonduplicate = copy.deepcopy(rows)
        for i, row in enumerate(nonduplicate):
            row["LOG2_FC_IP"] = row["CDS_GC"] + 1e-5 * ((i % 2) * 2 - 1)
        self.assertFalse(p.learn_transform(nonduplicate, "association")["association_projection"]["columns"][0]["zero"])

    def test_duplicate_model_cannot_change_ridge_penalty(self):
        rows = p.validate_features(self.train, "train")
        test = p.validate_features(self.holdout, "holdout")
        for row in rows + test:
            row["LOG2_FC_IP"] = row["CDS_GC"]
        y = p.ordered_targets(self.y, rows)
        baseline = p.fit_ridge_model(rows, y, test, "codon_core")
        duplicate = p.fit_ridge_model(rows, y, test, "association")
        self.assertEqual(baseline["selected_lambda"], duplicate["selected_lambda"])
        np.testing.assert_allclose(baseline["holdout_prediction"], duplicate["holdout_prediction"], atol=1e-10, rtol=1e-10)

    def test_baseline_nonlinearity_not_unique_association(self):
        rows = p.validate_features(self.train, "train")
        test = p.validate_features(self.holdout, "holdout")
        for row in rows + test:
            row["LOG2_FC_IP"] = math.log(row["HL_WT"]) ** 3
        settings = p.learn_transform(rows, "association")
        self.assertTrue(settings["association_projection"]["columns"][0]["zero"])
        y = np.asarray([math.log(r["HL_WT"]) ** 3 for r in rows])
        truth = np.asarray([math.log(r["HL_WT"]) ** 3 for r in test])
        fitted = p.fit_ridge_model(rows, y, test, "codon_core")
        pred = np.asarray(fitted["holdout_prediction"])
        self.assertLess(np.mean(np.abs(pred - truth)), .05)

    def test_ridge_matches_direct_normal_equations(self):
        x = np.asarray([[1., 0, 4], [2, 1, 4], [4, 3, 4], [7, 2, 4]])
        y = np.asarray([2., 0, 5, 3])
        prepared = p.ridge_prepare(x, y)
        fitted = p.ridge_parameters(prepared, 10.)
        z = (x - prepared["mean"]) / prepared["scale"]
        z[:, prepared["constant"]] = 0
        direct = np.linalg.solve(z.T @ z + 10 * np.eye(3), z.T @ (y - y.mean()))
        np.testing.assert_allclose(fitted["beta"], direct, atol=1e-14)
        np.testing.assert_allclose(p.ridge_predict(x, fitted), y.mean() + z @ direct, atol=1e-14)

    def test_null_and_zero_error_baseline_are_not_success(self):
        target = {r["GeneID"]: math.log(r["HL_WT"]) for r in self.holdout}
        null_training = {r["GeneID"]: math.log(r["HL_WT"]) for r in self.train}
        null_fit = p.predict(self.train, null_training, self.holdout, {"status": "synthetic_only"})
        result = p.evaluate(null_fit, target)
        self.assertFalse(result["primary_success"])
        row = result["primary_comparisons"]["no_change"]
        self.assertEqual(result["metrics"]["no_change"]["MAE"], 0.0)
        self.assertIsNone(row["relative_MAE_reduction"])
        self.assertFalse(row["meets_5_percent"])
        # Also force every frozen predictor to the exact null response.
        artificial = copy.deepcopy(self.result)
        truth = [target[g] for g in artificial["holdout_ids"]]
        for model in p.MODELS:
            artificial["models"][model]["holdout_prediction"] = truth
        result = p.evaluate(artificial, target)
        self.assertFalse(result["primary_success"])
        self.assertIsNone(result["influence"]["codon_relative_gain_range"])

    def test_shared_baseline_coupling_known_answer(self):
        latent = np.asarray([-1., -1, 1, 1])
        error = np.asarray([-1., 1, -1, 1])
        measured_baseline = latent + error
        measured_depletion = latent
        change = measured_depletion - measured_baseline
        slope = np.dot(measured_baseline, change) / np.dot(measured_baseline, measured_baseline)
        self.assertEqual(slope, -.5)
        pred_endpoint = .5 * measured_baseline
        np.testing.assert_array_equal(pred_endpoint - measured_depletion,
                                      (pred_endpoint - measured_baseline) - change)

    def test_fixed_ties_and_pooled_mae_not_fold_mean(self):
        self.assertEqual(p.choose_lambda({l: 1.0 for l in p.LAMBDAS}), 1000.)
        scores = {l: 2.0 for l in p.LAMBDAS}
        scores[.1], scores[100.] = 1.0, 1.0 + .5e-12
        self.assertEqual(p.choose_lambda(scores), 100.)
        model = self.result["models"]["association"]
        y = np.asarray(self.result["training_logKO"])
        for key, values in model["oof_predictions"].items():
            self.assertEqual(model["pooled_oof_mae"][key], float(np.mean(np.abs(np.asarray(values) - y))))

    def test_common_time_unit_rescaling_all_models(self):
        factor = 60.0
        train, test = copy.deepcopy(self.train), copy.deepcopy(self.holdout)
        for row in train + test:
            row["HL_WT"] *= factor
        target = {g: y + math.log(factor) for g, y in self.y.items()}
        shifted = p.predict(train, target, test, {"status": "synthetic_only"})
        for model in p.MODELS:
            a = np.asarray(self.result["models"][model]["holdout_prediction"])
            b = np.asarray(shifted["models"][model]["holdout_prediction"])
            np.testing.assert_allclose(b - math.log(factor), a, atol=2e-9, rtol=1e-9)
            if model in p.MODELS[2:]:
                self.assertEqual(self.result["models"][model]["selected_lambda"], shifted["models"][model]["selected_lambda"])
        evaluation = p.evaluate(shifted, {g: y + math.log(factor) for g, y in self.hy.items()})
        for model in p.MODELS:
            self.assertAlmostEqual(self.evaluation["metrics"][model]["MAE"], evaluation["metrics"][model]["MAE"], places=9)
        old = self.evaluation["kinetic"]["holdout"]["rows"]
        new = evaluation["kinetic"]["holdout"]["rows"]
        np.testing.assert_allclose([r["dimensionless_rate_contrast"] for r in old],
                                   [r["dimensionless_rate_contrast"] for r in new], atol=1e-9, rtol=1e-10)

    def test_group_hash_permutation_duplicates_and_missing(self):
        import hashlib
        group = "SYN00000_0"
        expected = int(hashlib.sha256(("cnot3-focus-v1|inner|" + group).encode()).hexdigest(), 16) % 5
        self.assertEqual(p.group_fold(group, "inner"), expected)
        self.assertEqual(p.validate_features(list(reversed(self.train)), "train"),
                         p.validate_features(self.train, "train"))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            p.validate_features(self.train + [self.train[0]], "bad")
        with self.assertRaisesRegex(ValueError, "leakage"):
            p.check_design(self.train, self.holdout + [self.train[0]], {"status": "synthetic_only"})
        bad = copy.deepcopy(self.holdout)
        bad[0]["group_ID"] = self.train[0]["group_ID"]
        with self.assertRaises(ValueError):
            p.check_design(self.train, bad, {"status": "synthetic_only"})
        with self.assertRaisesRegex(ValueError, "Target GeneIDs"):
            p.evaluate(self.result, dict(list(self.hy.items())[1:]))

    def test_real_coverage_gate_cannot_use_synthetic_minima(self):
        receipt = dict(status="root_authorized_prediction", reproduction_gate_passed=True,
                       reference_eligibility_passed=True, grouping_complete=True,
                       figure2_slam_overlap_count=len(self.train) + len(self.holdout),
                       eligible_count=len(self.train) + len(self.holdout))
        with self.assertRaisesRegex(ValueError, "Insufficient"):
            p.check_design(self.train, self.holdout, receipt)
        receipt["figure2_slam_overlap_count"] = 10000
        with self.assertRaisesRegex(ValueError, "coverage"):
            p.check_design(self.train, self.holdout, receipt)

    def test_metrics_ties_equal_groups_strata_and_signed_rows(self):
        metric = p.metrics(np.asarray([1., 3., 10.]), np.asarray([0., 1., 4.]), ["a", "a", "b"])
        self.assertEqual(metric["MAE"], 3.)
        self.assertEqual(metric["equal_component_MAE"], 3.75)
        self.assertEqual(p.spearman(np.array([1, 1, 2]), np.array([3, 3, 4])), 1.0)
        self.assertIsNone(p.spearman(np.array([1, 1]), np.array([2, 3])))
        self.assertEqual(np.searchsorted([1, 1, 2], [1, 1.5, 2, 3], side="left").tolist(), [0, 2, 2, 3])
        self.assertEqual(sum(x["n"] for x in self.evaluation["baseline_strata"]), len(self.holdout))
        for row in self.evaluation["rows"]:
            for error in row["models"].values():
                self.assertEqual(error["signed_error"], error["predicted"] - row["log_HL_KO"])

    def test_prediction_hash_and_commit_guard_precede_target_read(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            prediction = root / "predictions.json"
            artifact = copy.deepcopy(self.result)
            artifact["file_provenance"] = dict(implementation_sha256=p.sha256(p.__file__))
            p.write_json(prediction, artifact)
            args = argparse.Namespace(predictions=str(prediction), expected_predictions_sha256="0" * 64,
                                      commit_receipt=str(root / "missing_commit.json"),
                                      holdout_targets=str(root / "must_not_open_targets.json"), out=str(root / "out"))
            with mock.patch.object(p, "read_json", side_effect=AssertionError("No JSON should be read")):
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    p.run_evaluate(args)
            args.expected_predictions_sha256 = p.sha256(prediction)
            commit = root / "commit.json"
            p.write_json(commit, dict(status="synthetic_only", predictions_sha256="bad"))
            args.commit_receipt = str(commit)
            with self.assertRaisesRegex(ValueError, "Commit receipt"):
                p.run_evaluate(args)

    def test_duplicate_target_json_keys_and_nonfinite_json_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bad.json"
            path.write_text('{"gene": 1, "gene": 2}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
                p.read_json(path)
            path.write_text('{"gene": NaN}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Nonfinite JSON"):
                p.read_json(path)

    def test_fresh_directory_guard_and_full_synthetic_cli_stages(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inputs = dict(training_features=self.train, training_targets=self.y,
                          holdout_features=self.holdout, design_receipt={"status": "synthetic_only"})
            paths = {}
            for key, data in inputs.items():
                path = root / (key + ".json")
                p.write_json(path, data)
                paths[key] = str(path)
            protocol = root / "protocol.md"
            protocol.write_text("Synthetic test only. No source observations.", encoding="utf-8")
            paths["protocol"] = str(protocol)
            freeze = root / "freeze.json"
            p.write_json(freeze, dict(status="synthetic_only", implementation_sha256=p.sha256(p.__file__),
                                      inputs={k: p.sha256(v) for k, v in paths.items()}))
            args = argparse.Namespace(**paths, freeze=str(freeze), out=str(root / "predicted"))
            # Numerics have already run above; this verifies real serialization and gates.
            with mock.patch.object(p, "predict", return_value=copy.deepcopy(self.result)):
                p.run_predict(args)
            with self.assertRaises(FileExistsError):
                p.run_predict(args)
            target_path = root / "holdout_targets.json"
            p.write_json(target_path, self.hy)
            prediction = root / "predicted" / "predictions.json"
            commit = root / "commit.json"
            p.write_json(commit, dict(status="synthetic_only", predictions_sha256=p.sha256(prediction),
                                     holdout_targets_sha256=p.sha256(target_path)))
            score_args = argparse.Namespace(predictions=str(prediction),
                                             expected_predictions_sha256=p.sha256(prediction),
                                             commit_receipt=str(commit), holdout_targets=str(target_path),
                                             out=str(root / "scored"))
            p.run_evaluate(score_args)
            self.assertTrue(p.read_json(root / "scored" / "evaluation.json")["primary_success"])
            with self.assertRaises(FileExistsError):
                p.run_evaluate(score_args)


if __name__ == "__main__":
    unittest.main(verbosity=2)
