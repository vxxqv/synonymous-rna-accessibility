"""Refuse runtime/script drift before launching a frozen scientific phase."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

EXPECTED={'implementation':'CPython','python':'3.12.14','numpy':'2.5.4','Biopython':'1.88'}

def verify_runtime(actual):
    if actual!=EXPECTED:raise RuntimeError('Frozen initial-run runtime mismatch: '+json.dumps(actual,sort_keys=True))

def main():
    if not __debug__:raise RuntimeError('Optimized execution is refused')
    p=argparse.ArgumentParser()
    p.add_argument('--script',required=True);p.add_argument('--expected-sha256',required=True)
    p.add_argument('--receipt',required=True);p.add_argument('arguments',nargs=argparse.REMAINDER)
    a=p.parse_args();receipt=Path(a.receipt)
    if receipt.exists():raise FileExistsError('Guard receipt must be fresh')
    script=Path(a.script).resolve()
    digest=hashlib.sha256(script.read_bytes()).hexdigest()
    if digest!=a.expected_sha256:raise RuntimeError('Frozen phase script hash mismatch')
    import numpy
    import Bio
    actual={'implementation':platform.python_implementation(),'python':platform.python_version(),'numpy':numpy.__version__,'Biopython':Bio.__version__}
    verify_runtime(actual)
    arguments=a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
    receipt.parent.mkdir(parents=True,exist_ok=True)
    record={'status':'runtime_and_script_verified_before_phase','runtime':actual,'script_sha256':digest,'guard_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'phase_returncode':None}
    receipt.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
    result=subprocess.run([sys.executable,'-I','-B',str(script),*arguments])
    record['phase_returncode']=result.returncode
    record['status']='phase_completed' if result.returncode==0 else 'phase_failed_evidence_preserved'
    receipt.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
