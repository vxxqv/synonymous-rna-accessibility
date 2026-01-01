"""Write a complete synthetic-only prediction/evaluation receipt in a fresh folder."""
import argparse
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import predict as p
from test_predict import synthetic


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    train, y, holdout, hy = synthetic()
    inputs = dict(training_features=train, training_targets=y, holdout_features=holdout,
                  design_receipt={"status": "synthetic_only"})
    paths = {}
    for name, data in inputs.items():
        path = out / (name + ".json")
        p.write_json(path, data)
        paths[name] = str(path)
    protocol = out / "synthetic_protocol.md"
    protocol.write_text("Synthetic only. Deterministic seed741,180 repeated groups,2 members each. "
                        "Known response logKO=x+0.15*x^2-0.08*x^3+1.7*IP, x=logWT. "
                        "This is an implementation fixture, not biological evidence.\n", encoding="utf-8")
    paths["protocol"] = str(protocol)
    freeze = out / "synthetic_freeze.json"
    p.write_json(freeze, dict(status="synthetic_only", implementation_sha256=p.sha256(p.__file__),
                              inputs={key: p.sha256(value) for key, value in paths.items()}))
    p.run_predict(argparse.Namespace(**paths, freeze=str(freeze), out=str(out / "predictions")))
    target = out / "holdout_targets.json"
    p.write_json(target, hy)
    prediction = out / "predictions" / "predictions.json"
    commit = out / "synthetic_commit_receipt.json"
    p.write_json(commit, dict(status="synthetic_only", predictions_sha256=p.sha256(prediction),
                             holdout_targets_sha256=p.sha256(target)))
    p.run_evaluate(argparse.Namespace(predictions=str(prediction), expected_predictions_sha256=p.sha256(prediction),
                                      commit_receipt=str(commit), holdout_targets=str(target),
                                      out=str(out / "evaluation")))
    positive = p.read_json(out / "evaluation" / "evaluation.json")
    null_y = {r["GeneID"]: math.log(r["HL_WT"]) for r in train}
    null_hy = {r["GeneID"]: math.log(r["HL_WT"]) for r in holdout}
    null_prediction = p.predict(train, null_y, holdout, {"status": "synthetic_only"})
    null = p.evaluate(null_prediction, null_hy)
    p.write_json(out / "null_prediction.json", null_prediction)
    p.write_json(out / "null_evaluation.json", null)
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).parent), pattern="test_predict.py")
    with (out / "test_log.txt").open("x", encoding="utf-8") as log:
        tests = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    receipt = dict(status="synthetic_checks_passed" if tests.wasSuccessful() else "synthetic_checks_failed",
                   tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
                   source_sha256={name: p.sha256(Path(__file__).parent / name)
                                  for name in ("predict.py", "test_predict.py", "run_synthetic.py")},
                   python=platform.python_version(), numpy=np.__version__,
                   real_outcome_reads=0, real_fits=0, training_genes=len(train), holdout_genes=len(holdout),
                   positive_control_success=positive["primary_success"],
                   positive_control_MAE={model: item["MAE"] for model, item in positive["metrics"].items()},
                   positive_control_codon_gain=positive["primary_comparisons"]["codon_core"]["relative_MAE_reduction"],
                   null_control_success=null["primary_success"], null_no_change_MAE=null["metrics"]["no_change"]["MAE"],
                   null_relative_gain=null["primary_comparisons"]["no_change"]["relative_MAE_reduction"],
                   artifact_sha256={"predictions": p.sha256(prediction),
                                    "evaluation": p.sha256(out / "evaluation" / "evaluation.json"),
                                    "test_log": p.sha256(out / "test_log.txt")})
    p.write_json(out / "receipt.json", receipt)
    if not tests.wasSuccessful() or not positive["primary_success"] or null["primary_success"]:
        raise SystemExit("Synthetic falsification: inspect retained artifacts")
    print(receipt)


if __name__ == "__main__":
    main()
