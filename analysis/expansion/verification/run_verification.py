"""Portable final verifier with source hashes observed before and after execution.

Run from the repository root:
python analysis/expansion/verification/run_verification.py --outdir results/expansion/verification_final
This reruns the declared finite universes, never the large native folding study.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--outdir',required=True);parser.add_argument('--command-timeout',type=float,default=180);args=parser.parse_args()
    if not __debug__:raise RuntimeError('assertions must remain enabled')
    root=Path(__file__).resolve().parents[3];out=Path(args.outdir)
    if not out.is_absolute():out=root/out
    out=out.resolve()
    if root not in out.parents:raise ValueError('output must remain within repository')
    if out.exists():raise ValueError('fresh output directory required')
    out.mkdir(parents=True)
    paths=['analysis/expansion/verification/'+name for name in('run_verification.py','arbitrary_pair_graph_oracle.py','leucine_hardcore_volume.py','review_score_controls.py')]
    paths += ['analysis/exact/independent_oracle.py','analysis/expansion/graph_compiler.py','analysis/expansion/validate_graph_compiler.py','analysis/expansion/pilot_figures.py','analysis/expansion/conditional_vertex_null.py','analysis/expansion/retention_filtration.py','analysis/expansion/provenance.py','src/accessibility.py']
    before={p:digest(root/p)for p in paths};started=datetime.datetime.now(datetime.timezone.utc).isoformat();clock=time.perf_counter();commands=[];status='completed'
    specs=[('arbitrary_pair_graph_oracle.py',[sys.executable,str(root/'analysis/expansion/verification/arbitrary_pair_graph_oracle.py'),'--outdir',str(out)]),('leucine_hardcore_volume.py',[sys.executable,str(root/'analysis/expansion/verification/leucine_hardcore_volume.py'),'--outdir',str(out)]),('review_score_controls.py',[sys.executable,str(root/'analysis/expansion/verification/review_score_controls.py'),'--outdir',str(out)]),('validate_graph_compiler.py',[sys.executable,str(root/'analysis/expansion/validate_graph_compiler.py')]),('pilot_figures.py',[sys.executable,str(root/'analysis/expansion/pilot_figures.py')])]
    for name,argv in specs:
        timer=time.perf_counter()
        try:
            completed=subprocess.run(argv,cwd=root,capture_output=True,text=True,timeout=args.command_timeout)
            (out/(name+'.stdout.txt')).write_text(completed.stdout);(out/(name+'.stderr.txt')).write_text(completed.stderr)
            commands.append({'script':name,'argv':argv,'exit_status':completed.returncode,'seconds':time.perf_counter()-timer})
            if completed.returncode:status='failed';break
            print(name,'passed',round(commands[-1]['seconds'],3),'s',flush=True)
        except subprocess.TimeoutExpired as exc:
            status='timeout';commands.append({'script':name,'argv':argv,'exit_status':None,'seconds':time.perf_counter()-timer});break
    after={p:digest(root/p)for p in paths}
    changed=[p for p in paths if before[p]!=after[p]]
    if changed:status='inputs-changed'
    # API/figure scripts deliberately keep their normal reusable output locations;
    # copy their fresh tables into this immutable verification run as well.
    if status=='completed':
        for name in('graph_compiler_validation.json','general_graph_figure_sources.json'):(out/name).write_bytes((root/'results/expansion'/name).read_bytes())
    try:
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
        dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True))
    except (OSError,subprocess.CalledProcessError):commit='unavailable';dirty=True
    import numpy
    result_files=[p for p in out.iterdir()if p.is_file()]
    manifest={'schema_version':2,'claim_id':'GENERAL_PAIR_GRAPH_AND_ADDITIVE_CONTROL_PILOT','repository':{'commit':commit,'dirty':dirty},'command':[sys.executable,'analysis/expansion/verification/run_verification.py','--outdir',str(out.relative_to(root))],'environment':{'software':[{'name':'CPython','version':platform.python_version()},{'name':'NumPy','version':numpy.__version__}],'hardware':platform.platform()+' '+platform.machine()},'mathematics':{'assertion_tested':'Every nonempty1–2-codon/simple-H graph satisfies signature/geodesic/height3/order count assertions; all such signature fibres satisfy Leu-only activity2 volume law; finite control/barrier and API fixtures agree with independent exact checks.','coefficient_domain':'Z, Q and integer subset counts; ancillary float scores compared within1e-12','conventions':'standard61 sense codons, six ordered pair types, one-base moves, arbitrary simple simultaneous compatibility H; no physical multi-partner bonding','inputs':['20^m peptides and2^(choose(3m,2)) graphs for m1,2; all nonempty cases covered by H subset allowed-edges(genotype)','all fixed-base/partial-RY local codon cases','four-vertex graph/weight controls; complete capped API genotype spaces'],'bounds':{'all_graph_peptide_instances':13107360,'nonempty_instances_expected':162263,'connected_endpoint_pairs_expected':719338,'fixed_local_cases':70875,'additive_graph_weight_cases':5184},'non_claims':['Universal claims depend on the separate structural proofs.','No new native folding or biological-selection claim is established.','No bound on whole-ensemble folding kinetics or genotype fitness.']},'randomness':{'used':True,'generator':'PythonMT19937 and NumPy default_rng for finite test fixtures only','seed':[2026102402,2026102403,2026102404,42]},'run':{'started_at':started,'runtime_seconds':time.perf_counter()-clock,'exit_status':0 if status=='completed'else 1,'status':status,'resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':None,'max_threads':None},'per_command_timeout_seconds':args.command_timeout},'input_artifacts':[{'path':p,'sha256':before[p],'sha256_after':after[p]}for p in paths],'outputs':[{'path':str(p.relative_to(root)).replace('\\','/'),'sha256':digest(p),'kind':'result'if p.suffix=='.json'else'log'}for p in result_files],'commands':commands,'checks':[{'name':'input-hashes-unchanged','passed':not changed}],'result':status,'residual_risks':['No hard memory/CPU limits were imposed.','The independent code uses inherited model definitions and is not formal proof-assistant verification.'],'declared_results':[str(p.relative_to(root)).replace('\\','/')for p in result_files if p.suffix=='.json']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('final status:',status,'seconds:',round(manifest['run']['runtime_seconds'],3),flush=True)
    if status!='completed':sys.exit(1)

if __name__=='__main__':main()

