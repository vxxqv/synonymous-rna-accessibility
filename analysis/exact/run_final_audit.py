"""One bounded provenance run for the final independent exact computation suite."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);args=p.parse_args()
    if not __debug__:raise RuntimeError('assertions must be enabled')
    root=Path(__file__).resolve().parent
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    commands=['independent_oracle.py','register_fibonacci.py','bridge_transfer.py','peptide_automaton.py','bridge_threshold_automaton.py']
    observed=[]
    for name in commands:
        started=time.perf_counter()
        run=subprocess.run([sys.executable,str(root/name)],capture_output=True,text=True,check=False)
        (out/(name+'.stdout.txt')).write_text(run.stdout,encoding='utf-8')
        (out/(name+'.stderr.txt')).write_text(run.stderr,encoding='utf-8')
        observed.append({'script':name,'exit_code':run.returncode,'seconds':time.perf_counter()-started})
        if run.returncode:raise RuntimeError((name,run.returncode,run.stderr))
    names=['independent_results.json','register_examples.json','register_fibonacci_results.json','register_fibonacci_examples.json','transfer_validation.json','exact_peptide_ensemble.json','exact_peptide_automata.json','exact_bridge_threshold_automata.json']
    results={}
    for name in names:
        payload=(root/name).read_bytes();(out/name).write_bytes(payload)
        results[name]={'sha256':hashlib.sha256(payload).hexdigest(),'data':json.loads(payload)}
    report={'status':'all assertions passed','python':platform.python_version(),'assertions_enabled':__debug__,'commands':observed,'results':results,'nonclaims':['No thermodynamic energy or native biological result established.','Universal theorems require the accompanying proofs, not just these finite checks.','No random peptide-ensemble samples are run; exact automata supersede exploratory samples.'],'known_regressions_checked':{'empty_C0_nonempty_Cd':'minimum(K,K,3)=(None,0,0,4)','delta0_without_strict_mutation_bridge':'G(C,I,2)={0:1}'}}
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'seconds':sum(x['seconds'] for x in observed),'scripts':len(commands),'result_files':len(names)}))

if __name__=='__main__':main()
