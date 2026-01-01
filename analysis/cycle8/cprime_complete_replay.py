"""Replay the independent B1 source reconstruction from a portable archive."""
import argparse,subprocess,sys,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();root=a.repo.resolve();base=root/'analysis/cycle8'
subprocess.run([sys.executable,'-B',str(base/'cprime_independent_replay/replay_schema_adapter_v1.py'),
 '--source',str(base/'cprime_gate/supp_Supplementary Data 8.xlsx'),
 '--protocol',str(base/'cprime_pilot/protocol.md'),'--eligibility',str(base/'cprime_pilot/eligibility.json'),
 '--result',str(root/'results/cycle8/cprime/primary_result.json'),
 '--result-sha256','eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf',
 '--authorization',str(base/'cprime_independent_replay/actual_authorization.json'),
 '--output',str(a.out.resolve())],check=True)
