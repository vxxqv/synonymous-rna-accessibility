"""Run the sealed checker once; observe successful numeric-identity comparisons only.

The profiler never changes scientific functions, inputs, outputs or predicates.
"""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True


def main(authorization, out):
    source = Path(__file__).with_name('reference_binary_float.py')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != 'dd84524b3909df0437b070ebf1d97d7ba41449389dee78af53b0f1f86991e538':
        raise ValueError('Sealed reference changed')
    spec = importlib.util.spec_from_file_location('sealed_binary_float_reference', source)
    ref = importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
    observations = {}
    def profile(frame, event, arg):
        if event != 'return' or frame.f_code is not ref.assert_close.__code__:
            return
        a, b, label = [frame.f_locals[k] for k in ('a', 'b', 'location')]
        group = label.split('/')[0] if '/' in label else label
        entry = observations.setdefault(group, {'comparisons': 0, 'numeric_comparisons': 0,
            'maximum_absolute_difference': 0.0, 'maximum_scaled_difference': 0.0})
        entry['comparisons'] += 1
        if a is None or b is None or isinstance(a, (str, bool)) or isinstance(b, (str, bool)):
            return
        x, y = float(a), float(b)
        if not math.isfinite(x) or not math.isfinite(y):
            return
        error = abs(x-y)
        entry['numeric_comparisons'] += 1
        entry['maximum_absolute_difference'] = max(entry['maximum_absolute_difference'], error)
        entry['maximum_scaled_difference'] = max(entry['maximum_scaled_difference'], error/max(1.0,abs(x),abs(y)))
    error = None
    try:
        sys.setprofile(profile)
        result = ref.run(authorization, out)
    except Exception as exc:
        error = exc
    finally:
        sys.setprofile(None)
    destination = Path(out)
    if destination.exists():
        receipt = {'scope': 'Observation-only profiler of the unchanged sealed checker; no functions patched',
            'reference_sha256': digest, 'observer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'groups': observations, 'checker_completed': error is None}
        (destination/'numeric_identity_observations.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    if error is not None:
        raise error
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--authorization',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.authorization,a.out)
