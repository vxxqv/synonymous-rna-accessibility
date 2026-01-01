"""Small guard tests. No scientific outcomes are decoded or recomputed."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("archive_replay", HERE / "replay.py")
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class Guards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cnot3_replay_guards_")
        self.root = Path(self.temp.name)
        self.archive = self.root / "archive"
        self.archive.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def test_pinned_manifest(self):
        self.assertEqual(replay.digest(HERE / "manifest.json"), replay.MANIFEST_SHA256)

    def test_existing_output_refused_before_decode_or_child(self):
        output = self.root / "exists"
        output.mkdir()
        with patch.object(replay, "checked_inputs", side_effect=AssertionError("input access")), \
             patch.object(replay.subprocess, "run", side_effect=AssertionError("child launched")):
            with self.assertRaises(FileExistsError):
                replay.replay(self.archive, output)

    def test_output_inside_archive_refused_before_decode_or_child(self):
        with patch.object(replay, "checked_inputs", side_effect=AssertionError("input access")), \
             patch.object(replay.subprocess, "run", side_effect=AssertionError("child launched")):
            with self.assertRaisesRegex(ValueError, "outside the archive"):
                replay.replay(self.archive, self.archive / "new")

    def test_missing_file_refused_without_child(self):
        output = self.root / "missing_output"
        with patch.object(replay.subprocess, "run", side_effect=AssertionError("child launched")):
            with self.assertRaises(FileNotFoundError):
                replay.replay(self.archive, output)
        self.assertFalse(output.exists())

    def test_role_hash_drift_refused_before_scientific_decode(self):
        # Deliberately invalid JSON proves the guard works on bytes, not decoded
        # numerical input. The same routine checks both source and data roles.
        file = self.archive / "input.json"
        good = b"NOT SCIENTIFIC JSON"
        file.write_bytes(good)
        for role in ("training_features", "holdout_targets", "author_code", "verifier", "reference_core"):
            with self.subTest(role=role):
                manifest = {"inputs": {role: {"path": "input.json", "sha256": hashlib.sha256(good).hexdigest()}},
                            "supporting_evidence": {}}
                file.write_bytes(good)
                paths, hashes = replay.checked_inputs(self.archive, manifest)
                self.assertEqual(len(hashes), 1)
                self.assertEqual(paths[role].read_bytes(), good)
                file.write_bytes(good + b"!")
                with self.assertRaisesRegex(ValueError, "before scientific input decoding"):
                    replay.checked_inputs(self.archive, manifest)

    def test_supporting_evidence_hash_drift(self):
        file = self.archive / "review.md"
        file.write_text("review", encoding="utf-8")
        manifest = {"inputs": {}, "supporting_evidence": {
            "cohort_review": {"path": "review.md", "sha256": replay.digest(file)}}}
        file.write_text("changed", encoding="utf-8")
        with self.assertRaises(ValueError):
            replay.checked_inputs(self.archive, manifest)

    def test_path_escapes_refused(self):
        outside = self.root / "outside.txt"
        outside.write_text("outside", encoding="utf-8")
        for path in ("../outside.txt", "/outside.txt", "C:/outside.txt", "..\\outside.txt"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                replay.contained_file(self.archive, path)

    def test_duplicate_and_nonfinite_json_refused(self):
        path = self.root / "bad.json"
        for value in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}'):
            path.write_text(value, encoding="utf-8")
            with self.subTest(value=value), self.assertRaises(ValueError):
                replay.strict_json(path)

    def test_manifest_drift_refused_before_paths(self):
        staged = self.root / "wrapper"
        staged.mkdir()
        (staged / "replay.py").write_bytes((HERE / "replay.py").read_bytes())
        (staged / "manifest.json").write_bytes((HERE / "manifest.json").read_bytes() + b" ")
        result = subprocess.run([sys.executable, "-I", "-B", str(staged / "replay.py"),
                                 "--archive-root", str(self.root / "absent"), "--out", str(self.root / "new")],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("manifest SHA256 mismatch", result.stderr)
        self.assertFalse((self.root / "new").exists())

    def test_optimized_execution_refused_before_input_access(self):
        for flag in ("-O", "-OO"):
            result = subprocess.run([sys.executable, "-I", flag, str(HERE / "replay.py"),
                                     "--archive-root", str(self.root / "absent"), "--out", str(self.root / "new")],
                                    capture_output=True, text=True)
            with self.subTest(flag=flag):
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Optimized Python is refused before input access", result.stderr)
                self.assertFalse((self.root / "new").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
