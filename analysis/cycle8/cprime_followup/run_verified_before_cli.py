"""Safe entry point retaining the exact reviewed analysis source."""
if not __debug__:
    raise RuntimeError('Optimized Python is prohibited before any data access')
from pathlib import Path
import hashlib
import runpy

p=Path(__file__).with_name('analysis.py')
if hashlib.sha256(p.read_bytes()).hexdigest()!='45f75f86ef97325ae17bcabd306504dc57f5d885e81d3244a71eb43a5870619f':
    raise RuntimeError('The reviewed analysis source changed')
runpy.run_path(str(p),run_name='__main__')
