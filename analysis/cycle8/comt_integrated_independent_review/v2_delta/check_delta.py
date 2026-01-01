"""Narrow independent v1-to-v2 review, metadata and invented values only."""
import ast,difflib,hashlib,importlib.util,io,json,shutil,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Optimized execution refused')
P=Path(__file__).resolve().parent
BASE=P.parents[1]/'comt_pilot'
V1=BASE/'implementation_v1';V2=BASE/'implementation_v2'
H1='85f960e0ea4062cf013a655297c4c09108e7379ff133037c21f529b03e2a0cfd'
H2='f666dc9415d523b18350ffdedf5606e050803fa3e5c5da2efdd86f50f33e750d'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks=[]
def check(label,v):
    if not v:raise AssertionError(label)
    checks.append(label)
check('original manifest intact',sha(V1/'candidate_manifest.json')==H1)
check('successor manifest intact',sha(V2/'candidate_manifest_v2.json')==H2)
m1=json.loads((V1/'candidate_manifest.json').read_text());m2=json.loads((V2/'candidate_manifest_v2.json').read_text())
for folder,manifest in ((V1,m1),(V2,m2)):
    for n,h in manifest['files_sha256'].items():check(str(folder.name)+' pin '+n,sha(folder/n)==h)
    for n,h in manifest['held_sources_unchanged'].items():check(str(folder.name)+' held '+n,sha(BASE/n)==h)
check('same held protocol, folds, kernel and runtime',m1['held_sources_unchanged']==m2['held_sources_unchanged'])
changed=[p.name for p in V1.iterdir() if p.is_file() and (V2/p.name).exists() and p.read_bytes()!=(V2/p.name).read_bytes()]
check('only reader and test changed among original top level files',sorted(changed)==['selective_reader.py','test_implementation.py'])
added="            if n!='CAF' and v<0:raise ValueError('Negative required classifier feature '+n)\n"
r1=(V1/'selective_reader.py').read_text();r2=(V2/'selective_reader.py').read_text()
check('production delta exactly one domain guard',r2.count(added)==1 and r2.replace(added,'')==r1)
a1=ast.parse((V1/'test_implementation.py').read_text());a2=ast.parse((V2/'test_implementation.py').read_text())
for node in a2.body:
    if isinstance(node,ast.ClassDef):node.body=[n for n in node.body if not isinstance(n,ast.FunctionDef) or n.name!='test_classifier_features_nonnegative_domain']
check('test delta only the declared new test',ast.dump(a1)==ast.dump(a2))
diff=''.join(difflib.unified_diff(r1.splitlines(True),r2.splitlines(True),fromfile='v1reader',tofile='v2reader'))
(P/'independent_reader.diff').write_text(diff,encoding='utf8')
mirror=P/'mirror'
if mirror.exists():raise FileExistsError('Fresh mirror required')
mirror.mkdir()
for n in ('driver.py','selective_reader.py','normalize.py','prediction.py','test_implementation.py'):shutil.copyfile(V2/n,mirror/n)
shutil.copyfile(BASE/'proxy_models.py',P/'proxy_models.py')
sys.path.insert(0,str(mirror))
import selective_reader as new
s=importlib.util.spec_from_file_location('old_reader',V1/'selective_reader.py');old=importlib.util.module_from_spec(s);s.loader.exec_module(old)
from test_implementation import synthetic_row,small_contract,bytes_table
rows=[]
for name in new.NUMERIC[1:]:
    for value in (-1,0,1):
        row=synthetic_row();row[name]=value
        out={}
        for key,reader in [('old',old),('new',new)]:
            header,stream=bytes_table([row])
            try:
                parsed,_=reader.parse_summary(stream,header,{'ROI':1,'sample':'S'},small_contract())
                out[key]={'status':'accepted','value':parsed[0][name]}
            except ValueError as e:out[key]={'status':'rejected','error':str(e)}
        check('old accepts '+name+'='+str(value),out['old']['status']=='accepted')
        check('new correct domain '+name+'='+str(value),out['new']['status']==('rejected' if value<0 else 'accepted'))
        if value>=0:check('valid parsed value unchanged '+name+'='+str(value),out['old']==out['new'])
        rows.append({'feature':name,'input':value,**out})
f1=V1/'synthetic_stages_v2/evaluation/evaluation.json'
f2=V2/'synthetic_delta_stages_v1/evaluation/evaluation.json'
check('all125synthetic evaluation bytes identical',f1.read_bytes()==f2.read_bytes())
check('synthetic output matches declared checksum',sha(f2)=='30169a53447925632382361d0e0e9b9b2816b72a18110d3a77bd0f0c582f3ed9')
run=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s',str(mirror),'-p','test_implementation.py','-v'],capture_output=True,text=True)
(P/'unit_replay.log').write_text(run.stdout+'\n'+run.stderr,encoding='utf8')
check('all16unchangedcandidate unit tests pass',run.returncode==0 and 'Ran 16 tests' in run.stderr)
for folder,m in ((V1,m1),(V2,m2)):
    for n,h in m['files_sha256'].items():check('after '+folder.name+' '+n,sha(folder/n)==h)
receipt={'status':'PASS_NARROW_V2_DELTA','scope':'Closes v1 negative-feature guard HOLD only; root freeze/phase authorization still required. No real-source run or renewed biological validation.','predecessor_manifest_sha256':H1,'candidate_manifest_sha256':H2,'checks':checks,'check_count':len(checks),'domain_cases':rows,'unit_tests':16,'evaluation_sha256':sha(f2),'scientific_implementations_unchanged':['normalize.py','prediction.py','proxy_models.py','driver.py','apply_source_rf.R','source_helpers.R'],'real_source_rows_read':0,'source_simulation_rows_read':0,'python':sys.version,'verifier_sha256':sha(Path(__file__))}
(P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'units':16,'evaluation_sha256':sha(f2)},indent=2))
