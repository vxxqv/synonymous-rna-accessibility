"""Independent interface checks using invented numeric inputs and metadata only."""
import copy, hashlib, importlib.util, io, json, math, random, shutil, subprocess, sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

P=Path(__file__).resolve().parent
SRC=P.parent/'comt_pilot'/'implementation_v1'
PIN='85f960e0ea4062cf013a655297c4c09108e7379ff133037c21f529b03e2a0cfd'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')
if not __debug__:raise RuntimeError('Optimized Python refused before inputs')
if sha(SRC/'candidate_manifest.json')!=PIN:raise RuntimeError('Candidate changed')
manifest=json.loads((SRC/'candidate_manifest.json').read_text())
for name,h in manifest['files_sha256'].items():
    if sha(SRC/name)!=h:raise RuntimeError('Candidate file changed '+name)
for name,h in manifest['held_sources_unchanged'].items():
    if sha(SRC.parent/name)!=h:raise RuntimeError('Held source changed '+name)
mirror=P/'mirror'
if mirror.exists():raise FileExistsError('Fresh independent mirror required')
mirror.mkdir()
for name in ('driver.py','selective_reader.py','normalize.py','prediction.py','test_implementation.py'):
    shutil.copyfile(SRC/name,mirror/name)
shutil.copyfile(SRC.parent/'proxy_models.py',P/'proxy_models.py')
sys.path.insert(0,str(mirror))
import selective_reader as sr, normalize as nm, prediction as pr, driver
from test_implementation import synthetic_row,small_contract,bytes_table
checks=[];bugs=[]
def check(label,condition,evidence=None):
    if not condition:raise AssertionError(label)
    checks.append({'check':label,'passed':True,'evidence':evidence})
def rejected(fn):
    try:fn()
    except (ValueError,FileExistsError,RuntimeError):return True
    return False

# Source-only builder must regenerate all biological and target metadata.
s=importlib.util.spec_from_file_location('metadata_builder',SRC/'build_contract.py')
b=importlib.util.module_from_spec(s);s.loader.exec_module(b)
contract=json.loads((SRC/'selective_contract.json').read_text())
check('metadata contract independently rebuilt from pinned source and metadata',json.loads(json.dumps(b.build()))==contract)
check('63 required biological members',len(contract['members'])==63)
check('125 targets at53sites and113SNVs at52sites',len(contract['targets'])==125 and len({r['site'] for r in contract['targets']})==53 and sum(r['is_SNV'] for r in contract['targets'])==113 and len({r['site'] for r in contract['targets'] if r['is_SNV']})==52)
check('held codon alternatives share fold',all(len({r['assessment_fold'] for r in contract['targets'] if r['site']==site})==1 for site in {r['site'] for r in contract['targets']}))

# Independent exact normalization oracle: random two-base duplicate rows and
# a nontarget allele must enter the same WT denominator after median collapse.
rng=random.Random(84819)
for i in range(50):
    count=rng.randrange(1,30);other=rng.randrange(1,30);nc=rng.randrange(0,3)
    dp=200;rows=[synthetic_row(ref='AA',alt='CC',dp=dp-2,cao=count-1,posnt=1),synthetic_row(ref='AA',alt='CC',dp=dp+2,cao=count+1,posnt=2,line=3),synthetic_row(alt='G',dp=dp,cao=other,line=4),synthetic_row(sample='NC',ref='AA',alt='CC',dp=dp,cao=nc)]
    c=small_contract();c['targets'][0]['edit_key']='1:A>C;2:A>C'
    st={r['record_id']:{'VAR_ID':nm.raw_id(r),'frequency_pass':True,'globally_rejected':False} for r in rows}
    expected=F(2*count+1,2*(dp-count-other)+1)-F(2*nc+1,2*(dp-nc)+1)
    m=nm.normalized(rows,st,c)['measurements'][0]
    check('exact denominator andNCcase'+str(i),F(m['corrected_q'])==expected and (m['log_q']==math.log(float(expected)) if expected>0 else m['log_q'] is None))

# Poisoned excluded data, domain failures, and explicit guard counterexamples.
for n in sr.NUMERIC[1:]:
    row=synthetic_row();row[n]=-1;header,stream=bytes_table([row])
    if not rejected(lambda:sr.parse_summary(stream,header,{'ROI':1,'sample':'S'},small_contract())):
        bugs.append({'id':'negative_feature_domain','field':n,'value':-1,'observed':'accepted','required':'stop before RF per integration addendum'})
for n in sr.NUMERIC:
    for v in (float('nan'),float('inf')):
        row=synthetic_row();row[n]=v;header,stream=bytes_table([row])
        check('nonfinite rejected '+n+str(v),rejected(lambda:sr.parse_summary(stream,header,{'ROI':1,'sample':'S'},small_contract())))

# Synthetic normalization field and missingness changes. Median of all logs is
# not median of technical-group medians when technical multiplicities differ.
c={'members':{},'targets':[{'hgvs_nt':'v','site':'1','ROI':1,'assessment_fold':0,'is_SNV':True}]};ms=[]
for trait in ('D','R','F','H','M'):
    for bio in range(1,5):
        for tech in range(3 if bio==1 else 1):
            sid=f'{trait}{bio}_{tech}';c['members'][sid]={'sample':sid,'ROI':1,'trait':trait,'code_bio':bio,'geo_bio':bio}
            val=[0,100,200][tech] if bio==1 else bio
            ms.append({'hgvs_nt':'v','Sample':sid,'log_q':val+(10 if trait=='R' else 20 if trait=='H' else 0)})
ss=nm.contrast_summaries({'measurements':ms},c)
check('technical then group medians independently weighted',ss['code'][0]['assays']['D']['center']==3.5)
check('sign and shared contrast algebra',ss['code'][0]['contrasts']['R']['center']==10 and ss['code'][0]['contrasts']['L']['center']==10 and ss['code'][0]['contrasts']['T']['center']==0)
for r in ms:
    if r['Sample']=='H4_0':r['log_q']=None
ss=nm.contrast_summaries({'measurements':ms},c)
check('missingload leaves RNA and target complete',ss['code'][0]['contrasts']['R']['complete'] and ss['code'][0]['contrasts']['T']['complete'] and not ss['code'][0]['contrasts']['L']['complete'])

# Evaluate fully invented sealed predictions against a separately calculated
# site-first MAE. No numerical fitting oracle is used for expected values.
k=pr.load_kernel(P/'proxy_models.py')
rows=[];pred=[]
for i,t in enumerate(contract['targets']):
    rows.append({**t,'contrasts':{n:{'complete':True,'center':0.,'resolved_sign':0} for n in ('R','T','L','L_mid')}})
    e=F(1+(i%9),10)
    pred.append({**{n:t[n] for n in ('hgvs_nt','site','ROI')},'fold':t['assessment_fold'],'predictions':{n:{'exact_stored_input_arithmetic':str(e if n=='RNA_load_linear' else 1)} for n in pr.NAMES}})
summaries={'code':rows,'geo':copy.deepcopy(rows)}
predictions={case:{mode:{'predictions':copy.deepcopy(pred),'unavailable':[]} for mode in ('heavy','mid')} for case in ('code','geo')}
ev=pr.evaluate(predictions,summaries,k)
groups={s:[p for p in pred if p['site']==s] for s in {p['site'] for p in pred}}
oracle=sum(sum(F(p['predictions']['RNA_load_linear']['exact_stored_input_arithmetic']) for p in pp)/len(pp) for pp in groups.values())/len(groups)
check('sitebalanced error exact independent mean',F(ev['cases']['code']['heavy']['scores']['site_MAE']['RNA_load_linear']['exact_stored_input_arithmetic'])==oracle,str(oracle))
snv={r['hgvs_nt'] for r in contract['targets'] if r['is_SNV']}
check('SNVsubset uses retained same predictions',ev['cases']['code']['heavy']['sensitivities']['SNV']['variants']==len(snv))
check('bothmaps known signal passes',ev['primary_operational_expansion_pass'])
for p in predictions['geo']['heavy']['predictions']:p['predictions']['RNA_load_linear']['exact_stored_input_arithmetic']='91/100'
check('one map below10percent fails primary',not pr.evaluate(predictions,summaries,k)['primary_operational_expansion_pass'])
for p in predictions['geo']['heavy']['predictions']:p['predictions']['RNA_load_linear']['exact_stored_input_arithmetic']='9/10'
check('exact10percent boundary passes',pr.evaluate(predictions,summaries,k)['primary_operational_expansion_pass'])
for case in predictions:
    for p in predictions[case]['heavy']['predictions']:p['predictions']['ROI_median']['exact_stored_input_arithmetic']='0'
check('zero baseline fails despite other positive gains',not pr.evaluate(predictions,summaries,k)['primary_operational_expansion_pass'])

run=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s',str(mirror),'-p','test_implementation.py','-v'],capture_output=True,text=True)
(P/'author_unit_replay.log').write_text(run.stdout+'\n'+run.stderr,encoding='utf8')
check('unchanged author15synthetic tests replayed inownedmirror',run.returncode==0)
for name,h in manifest['files_sha256'].items():
    if sha(SRC/name)!=h:raise RuntimeError('Candidate changed during review '+name)
receipt={'status':'HOLD_REQUIRED_GUARD_CORRECTION','assertion_tested':'Declared metadata selection, exact normalization and site-scoring behavior on invented inputs only','candidate_manifest_sha256':PIN,'independent_checks':len(checks),'checks':checks,'counterexamples':bugs,'source_rows_decoded':0,'source_simulation_rows_decoded':0,'author_tests':15,'python':sys.version,'reviewer_source_sha256':sha(__file__),'candidate_before_after_unchanged':True,'limits':'No real-source run, external replication, runtime RF rerun or biological validation. Previously reviewed kernel not rederived.'}
dump(P/'receipt.json',receipt)
print(json.dumps({'status':receipt['status'],'independent_checks':len(checks),'counterexamples':bugs},indent=2))
