"""Only synthetic diagnostics/configuration objects. No folding or observed rows."""
import sys
if not __debug__:
    raise RuntimeError('Disabled assertions are refused before input access')
sys.dont_write_bytecode = True
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import types
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent

def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

observer = load('observe_failure')
pre = load('config_preflight')

class SyntheticChecks(unittest.TestCase):
    def test_collector_full_streams_and_whitelist(self):
        with tempfile.TemporaryDirectory(prefix='persist-diagnostic-') as folder:
            root = Path(folder)
            (root / 'STARTED.json').write_text('{"rows":37}', encoding='utf-8')
            worker = root / 'synthetic_worker'; worker.mkdir()
            raw = b'Traceback\n::error::inert\n\xff\n'
            (worker / 'stderr.txt').write_bytes(raw)
            (worker / 'stdout.txt').write_bytes(b'')
            (worker / 'row.json').write_text('{"python_warnings":["synthetic"]}', encoding='utf-8')
            (worker / 'Turner2004.exported.par').write_bytes(b'synthetic')
            (worker / 'unrelated_private.txt').write_text('DO_NOT_READ', encoding='utf-8')
            result = observer.observe(root)
            entries = {item['file']: item for item in result['files']}
            self.assertEqual(len(entries), 5)
            self.assertEqual(entries['synthetic_worker/stderr.txt']['sha256'], hashlib.sha256(raw).hexdigest())
            self.assertEqual(entries['synthetic_worker/stderr.txt']['text_utf8_with_replacement'], raw.decode('utf-8', 'replace'))
            self.assertEqual(entries['synthetic_worker/stdout.txt']['bytes'], 0)
            self.assertNotIn('text_utf8_with_replacement', entries['synthetic_worker/Turner2004.exported.par'])
            cli = subprocess.run([sys.executable, '-I', '-B', str(HERE / 'observe_failure.py'), '--root', str(root)], capture_output=True, text=True)
            self.assertEqual(cli.returncode, 0)
            self.assertEqual(len(cli.stdout.splitlines()), 1)
            self.assertEqual(json.loads(cli.stdout), result)
            self.assertTrue(cli.stdout.startswith('{'))

    def test_missing_root_is_explicit(self):
        self.assertEqual(observer.observe(HERE / 'does_not_exist_synthetic')['status'], 'diagnostic_directory_absent')

    def fixture(self, root, body='return RNA, object()', parameter=b'synthetic', snapshot='{}'):
        repo = root / 'repo'; repo.mkdir()
        path = repo / pre.PREFIX; path.mkdir(parents=True)
        author = ('import RNA\nimport warnings\n'
                  'def configured_md(settings):\n    ' + body + '\n'
                  'def snapshot_md(md, RNA):\n    return ' + snapshot + '\n'
                  'def fold_one(*args):\n    raise AssertionError("FOLD_FORBIDDEN")\n'
                  'def sequence_features(*args):\n    raise AssertionError("FEATURES_FORBIDDEN")\n')
        (path / 'persist_prediction.py').write_text(author, encoding='utf-8')
        (path / 'protocol.json').write_text('{"vienna":{"effective_md":{}}}', encoding='utf-8')
        (path / 'Turner2004.par').write_bytes(parameter)
        pins = {name: pre.digest(repo / name) for name in pre.PINNED}
        fake_rna = types.SimpleNamespace(__version__='2.7.2', pf_float_precision=lambda: 0)
        def export(filename):
            Path(filename).write_bytes(b'synthetic')
            return 1
        fake_rna.params_save = export
        return repo, pins, fake_rna

    def run_preflight_case(self, body='return RNA, object()', parameter=b'synthetic', snapshot='{}', error=None, expected_stage='complete_without_folding'):
        with tempfile.TemporaryDirectory(prefix='persist-diagnostic-') as folder:
            root = Path(folder); repo, pins, rna = self.fixture(root, body, parameter, snapshot)
            out = root / 'out'
            with patch.object(pre, 'PINNED', pins), patch.dict(sys.modules, {'RNA': rna, 'numpy': types.SimpleNamespace(__version__='2.5.4')}):
                if error:
                    with self.assertRaises(error): pre.preflight(repo, out)
                else:
                    pre.preflight(repo, out)
                    with self.assertRaises(ValueError): pre.preflight(repo, out)
            report = json.loads((out / 'preflight.json').read_text(encoding='utf-8'))
            self.assertEqual(report['stage'], expected_stage)
            self.assertEqual(report['folds'], 0)
            self.assertEqual(report['sequence_rows_read'], 0)
            self.assertTrue(report['inputs_immutable'])
            return report

    def test_configuration_success_no_fold(self):
        self.assertEqual(self.run_preflight_case()['status'], 'configuration_only_passed')

    def test_original_exception_exposed(self):
        report = self.run_preflight_case(body='raise KeyError("synthetic actual exception")', error=KeyError, expected_stage='configured_md')
        self.assertIn('synthetic actual exception', report['traceback'])
        self.assertEqual(report['status'], 'failed')

    def test_warning_not_ignored(self):
        self.run_preflight_case(body='warnings.warn("synthetic warning"); return RNA, object()', error=UserWarning, expected_stage='configured_md')

    def test_snapshot_and_parameter_mismatch(self):
        self.run_preflight_case(snapshot='{"different":1}', error=RuntimeError, expected_stage='snapshot_md')
        self.run_preflight_case(parameter=b'different frozen bytes', error=RuntimeError, expected_stage='parameter_export')

    def test_optimization_refused_before_any_io(self):
        for script in ('observe_failure.py', 'config_preflight.py'):
            for flag in ('-O', '-OO'):
                proc = subprocess.run([sys.executable, flag, '-I', '-B', str(HERE / script)], capture_output=True, text=True)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn('Disabled assertions are refused before input access', proc.stderr)

if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SyntheticChecks))
    receipt = {'scope': 'synthetic only; no real configuration or folding was executed',
               'tests': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
               'status': 'passed' if result.wasSuccessful() else 'failed',
               'folds': 0, 'outcome_reads': 0}
    with (HERE / 'synthetic_receipt.json').open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2); stream.write('\n')
    sys.exit(0 if result.wasSuccessful() else 1)
