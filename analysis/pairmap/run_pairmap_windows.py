from pathlib import Path
import sys,json,time,math,collections
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];root=ROOT/'data/pairmap';sys.path.insert(0,str(ROOT/'analysis'));import native_models

def run():
 blocks=json.loads((root/'eligible_windows.json').read_text());outdir=ROOT/'results/pairmap';outdir.mkdir(exist_ok=True);tasks=[(b,c) for b in blocks for c in ['incell','cellfree'] if not (outdir/f"{b['id']}_{c}.json").exists()];print('pending',len(tasks),flush=True)
 with ProcessPoolExecutor(max_workers=2) as p:
  fs={p.submit(native_models.evaluate_one,t):t for t in tasks}
  for k,f in enumerate(as_completed(fs),1):
   r=f.result();(outdir/f"{r['id']}_{r['condition']}.json").write_text(json.dumps(r,allow_nan=False));print('completed',k,len(tasks),r['id'],r['condition'],r['n_compatible'],r['seconds'],flush=True)
if __name__=='__main__':run()
