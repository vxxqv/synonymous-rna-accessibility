"""Content-addressed validity checks for exploratory result caches."""
from pathlib import Path
from functools import lru_cache
import json,hashlib,platform
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

@lru_cache(maxsize=None)
def code_manifest(method):
    paths=[Path(__file__),Path(__file__).parent/method,ROOT/'src/accessibility.py']
    for name in ('retention_filtration.py','conditional_vertex_null.py','smoothness_controls.py'):
        path=Path(__file__).parent/name
        if path not in paths:paths.append(path)
    return {str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in paths}

def signature(method,block,condition,parameters):
    manifest={'method':method,'parameters':parameters,'source_result_sha256':digest(ROOT/'results/native'/f"{block['id']}_{condition}.json"),
      'block_sha256':hashlib.sha256(json.dumps(block,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
      'code_sha256':code_manifest(method),'python':platform.python_version(),'numpy':np.__version__}
    manifest['signature']=hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return manifest

def cache_valid(path,manifest):
    path=Path(path)
    if not path.exists():return False
    try:return json.loads(path.read_text()).get('run_provenance',{}).get('signature')==manifest['signature']
    except (ValueError,OSError):return False
