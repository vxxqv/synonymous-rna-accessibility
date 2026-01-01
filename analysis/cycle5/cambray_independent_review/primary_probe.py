"""Isolated synthetic invocation of the held primary API, never the reference."""
import sys
if not __debug__:raise RuntimeError('Disabled assertions refused')
sys.dont_write_bytecode=True
import hashlib
import importlib.util
import json
from pathlib import Path

path=Path(sys.argv[1]);expected=sys.argv[2]
if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Primary review point changed')
spec=importlib.util.spec_from_file_location('held_primary',path)
primary=importlib.util.module_from_spec(spec);spec.loader.exec_module(primary)
out=[]
for case in json.load(sys.stdin):
    operation=case['operation']
    if operation=='classify':result=primary.classify(case['graph'],case['values'],case['endpoints'],case.get('thresholds'),case.get('tolerance',0))
    elif operation=='interval':result=primary.interval_certificate(case['graph'],case['values'],case['endpoints'],case['thresholds'])
    elif operation=='projection':result=primary.additive_projection(case['values'])
    elif operation=='prediction':result=primary.missing_corner_predictions({int(k):v for k,v in case['known'].items()},case['target'])
    elif operation=='geometry':result=primary.threshold_geometry(case['graph'],case['values'],case['endpoints'])
    elif operation=='preparations':
        result=[primary.classify(case['graph'],[[row[i],row[4+j]] for row in case['values']],case['endpoints']) for i in range(4) for j in range(3)]
    else:raise ValueError(operation)
    out.append(result)
json.dump(out,sys.stdout,default=str,allow_nan=False)
