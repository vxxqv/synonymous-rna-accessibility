"""Frozen exploratory local-descriptor order tests with exact certificates."""
from __future__ import annotations
import argparse, collections, hashlib, json, platform, time
from fractions import Fraction as F
from pathlib import Path

MODELS=['strict_shared','context_specific_additive','context_specific_all_observed_pairs']
CELLS=['K562','KBM7']; EXONS=[17,25,27,39,42]; GAMMAS=[F(1),F(5,4),F(3,2),F(2)]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,default=str)+'\n',encoding='utf-8')

def skip(s,i):
    while i<len(s) and s[i].isspace():i+=1
    if i==len(s):raise ValueError('Truncated JSON')
    if s[i]=='"':
        i+=1
        while i<len(s):
            if s[i]=='\\':i+=2
            elif s[i]=='"':return i+1
            else:i+=1
        raise ValueError('Unterminated JSON string')
    if s[i] in '[{':
        close=']' if s[i]=='[' else '}';i+=1
        while i<len(s):
            if s[i].isspace() or s[i] in ',:':i+=1
            elif s[i]==close:return i+1
            else:i=skip(s,i)
        raise ValueError('Unterminated JSON container')
    j=i
    while j<len(s) and s[j] not in ',]} \t\r\n':j+=1
    if j==i:raise ValueError('Invalid JSON token')
    return j

def members(s):
    s=s.strip()
    if not s.startswith('{') or not s.endswith('}'):raise ValueError('Expected JSON object')
    i=1;seen=set();out={}
    while i<len(s)-1:
        while s[i].isspace() or s[i]==',':i+=1
        if s[i]=='}':break
        j=skip(s,i);key=json.loads(s[i:j]);i=j
        if not isinstance(key,str) or key in seen:raise ValueError('Duplicate or invalid JSON key')
        seen.add(key)
        while s[i].isspace():i+=1
        if s[i]!=':':raise ValueError('Expected colon')
        i+=1
        while s[i].isspace():i+=1
        j=skip(s,i);out[key]=s[i:j];i=j
    return out

def elements(s):
    s=s.strip()
    if not s.startswith('[') or not s.endswith(']'):raise ValueError('Expected JSON array')
    i=1
    while i<len(s)-1:
        while s[i].isspace() or s[i]==',':i+=1
        if s[i]==']':break
        j=skip(s,i);yield s[i:j];i=j

def selected(s,fields):
    obj=members(s)
    if not set(fields)<=set(obj):raise ValueError('Missing selected source fields')
    return {k:json.loads(obj[k]) for k in fields}

def parse_b1(path):
    top=members(path.read_text(encoding='utf-8'))
    status=json.loads(top['status'])
    if status!='COMPLETE_DESCRIPTIVE_AUDIT':raise ValueError('B1 source audit did not complete')
    state_fields=['state','cell','usable_rows','complete','lower','upper','source_rows']
    states=[selected(s,state_fields) for s in elements(top['state_cells'])]
    source=[selected(s,['source_row','state','cell','exon']) for s in elements(top['source_rows'])]
    sensitivities=[]
    for raw in elements(top['sensitivities']):
        obj=members(raw)
        sensitivities.append({'gamma':json.loads(obj['gamma']),
          'exon_shared_order_results':[selected(s,['exon','status']) for s in elements(obj['exon_shared_order_results'])]})
    return states,source,sensitivities

def exact_solve(rows,rhs):
    """A rational solution with free variables zero, or refusal if inconsistent."""
    if not rows or not rows[0]:raise ValueError('Empty exact system')
    n=len(rows[0]);a=[[F(v) for v in row]+[F(b)] for row,b in zip(rows,rhs)];k=0;piv=[]
    for j in range(n):
        i=next((i for i in range(k,len(a)) if a[i][j]),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k];z=a[k][j];a[k]=[v/z for v in a[k]]
        for t in range(len(a)):
            if t!=k and a[t][j]:
                z=a[t][j];a[t]=[x-z*y for x,y in zip(a[t],a[k])]
        piv.append(j);k+=1
        if k==len(a):break
    if any(not any(row[:n]) and row[n] for row in a):raise ArithmeticError('Inconsistent rational support')
    ans=[F(0)]*n
    for i,j in enumerate(piv):ans[j]=a[i][n]
    if any(sum(x*y for x,y in zip(row,ans))!=b for row,b in zip(rows,rhs)):raise ArithmeticError('Exact solve mismatch')
    return ans

def validate_design(X,records,intervals):
    if not X or not X[0] or len(X)!=len(records) or len(X)!=len(intervals):raise ValueError('Nonempty equal-length inputs required')
    if len({len(r) for r in X})!=1 or any(type(v) is not int for row in X for v in row):raise ValueError('Integral rectangular matrix required')
    if len({(r['state'],r['cell'],r['exon']) for r in records})!=len(records):raise ValueError('Duplicate analysis row')
    for lo,up in intervals:
        if not isinstance(lo,F) or (up is not None and not isinstance(up,F)):raise ValueError('Rational intervals required')
        if lo<0 or (up is not None and (up<=0 or up<lo)):raise ValueError('Invalid positive-domain interval')

def forced(X,records,intervals):
    validate_design(X,records,intervals);pairs=[];A=[]
    for i,(li,ui) in enumerate(intervals):
        for j,(lj,uj) in enumerate(intervals):
            if records[i]['cell']==records[j]['cell'] and records[i]['exon']==records[j]['exon'] and uj is not None and li>uj:
                pairs.append((i,j));A.append([v-w for v,w in zip(X[i],X[j])])
    return pairs,A

def responses(records,intervals,latent):
    groups=collections.defaultdict(list)
    for i,r in enumerate(records):groups[r['cell'],r['exon']].append(i)
    values=[None]*len(records)
    for key,inds in groups.items():
        positives=[v for i in inds for v in intervals[i] if v is not None and v>0]
        y=min(positives)/2 if positives else F(1)
        for level in sorted({latent[i] for i in inds}):
            equal=[i for i in inds if latent[i]==level]
            y=max(y,max(intervals[i][0] for i in equal))
            if y<=0 or any(intervals[i][1] is not None and y>intervals[i][1] for i in equal):raise ArithmeticError('Positive monotone completion failed')
            for i in equal:values[i]=y
    return values

def verify_certificate(X,records,intervals,result):
    pairs,A=forced(X,records,intervals)
    if result['forced_pairs']!=[list(p) for p in pairs]:raise ArithmeticError('Certificate pair coverage mismatch')
    if result['status']=='feasible':
        beta=list(map(F,result['beta']));latent=list(map(F,result['latent']));y=list(map(F,result['responses']))
        if len(beta)!=len(X[0]) or len(latent)!=len(X) or len(y)!=len(X):raise ArithmeticError('Witness length')
        if latent!=[sum(v*b for v,b in zip(row,beta)) for row in X]:raise ArithmeticError('Latent mismatch')
        margins=[sum(v*b for v,b in zip(row,beta)) for row in A]
        if margins!=list(map(F,result['margins'])) or any(v<1 for v in margins):raise ArithmeticError('Margin failure')
        for i,(lo,up) in enumerate(intervals):
            if y[i]<=0 or y[i]<lo or (up is not None and y[i]>up):raise ArithmeticError('Interval response failure')
            for j in range(len(X)):
                if records[i]['cell']==records[j]['cell'] and records[i]['exon']==records[j]['exon'] and latent[i]<=latent[j] and y[i]>y[j]:raise ArithmeticError('Monotonicity failure')
    elif result['status']=='infeasible':
        weights=[F(0)]*len(A);seen=set()
        for term in result['dual']:
            k=term['constraint_index'];w=F(term['weight'])
            if type(k) is not int or not 0<=k<len(A) or k in seen or w<=0 or term['difference']!=A[k]:raise ArithmeticError('Dual support failure')
            seen.add(k);weights[k]=w
        if sum(weights)!=1 or any(sum(w*row[j] for w,row in zip(weights,A)) for j in range(len(X[0]))):raise ArithmeticError('Exact dual cancellation failure')
    else:raise ArithmeticError('Uncertified status')
    return True

def solve(X,records,intervals):
    """Floating candidate search, with exact rational evidence as sole acceptance."""
    import numpy as np
    from scipy.optimize import linprog
    pairs,A=forced(X,records,intervals);result={'forced_pairs':[list(p) for p in pairs],'forced_pair_count':len(pairs)}
    if not A:beta=[F(0)]*len(X[0]);optstatus='no_forced_pairs'
    else:
        opt=linprog(np.zeros(len(X[0])),A_ub=-np.array(A,dtype=float),b_ub=-np.ones(len(A)),bounds=[(None,None)]*len(X[0]),method='highs')
        optstatus=int(opt.status)
        if opt.status==0:
            beta=[F(float(v)).limit_denominator(10**6) for v in opt.x]
            gaps=[sum(v*b for v,b in zip(row,beta)) for row in A]
            if min(gaps)<=0:raise ArithmeticError('Candidate has no exact strict certificate')
            beta=[b/min(gaps) for b in beta]
        elif opt.status==2:
            ae=np.vstack([np.array(A,dtype=float).T,np.ones(len(A))]);rhs=np.r_[np.zeros(len(X[0])),1.]
            dual=linprog(np.zeros(len(A)),A_eq=ae,b_eq=rhs,bounds=(0,None),method='highs')
            if dual.status!=0:raise ArithmeticError('No dual candidate')
            support=[i for i,v in enumerate(dual.x) if v>1e-9]
            mat=[[A[i][j] for i in support] for j in range(len(X[0]))]+[[1]*len(support)]
            ws=exact_solve(mat,[F(0)]*len(X[0])+[F(1)])
            if any(w<0 for w in ws):raise ArithmeticError('Rational dual is not nonnegative')
            result.update(status='infeasible',dual=[{'constraint_index':i,'weight':str(w),'difference':A[i]} for i,w in zip(support,ws) if w],optimizer_candidate_status=optstatus)
            verify_certificate(X,records,intervals,result);result['exact_verified']=True;return result
        else:raise ArithmeticError('Unresolved optimizer status '+str(opt.status))
    latent=[sum(v*b for v,b in zip(row,beta)) for row in X]
    result.update(status='feasible',beta=list(map(str,beta)),latent=list(map(str,latent)),
      responses=list(map(str,responses(records,intervals,latent))),
      margins=[str(sum(v*b for v,b in zip(row,beta))) for row in A],optimizer_candidate_status=optstatus)
    verify_certificate(X,records,intervals,result);result['exact_verified']=True;return result

def load_contract(source,matrix_path,mapping_path):
    states,source_rows,sensitivities=parse_b1(source)
    matrices=json.loads(matrix_path.read_text());mapping=json.loads(mapping_path.read_text());metadata={m['id']:m for m in mapping}
    if len(metadata)!=66 or sorted({m['exon'] for m in mapping})!=EXONS:raise ValueError('Fixed metadata universe changed')
    expected=[{'cell':cell,'id':key,'exon':metadata[key]['exon']} for cell in CELLS for key in sorted(metadata)]
    if matrices['rows']!=expected or set(matrices['matrices'])!=set(MODELS):raise ValueError('Frozen matrix row or model mismatch')
    interval_by={};source_index={};by=collections.defaultdict(list)
    for row in source_rows:
        key=(row['state'],row['cell'])
        if row['state'] not in metadata or row['cell'] not in CELLS or type(row['exon']) is not int or row['exon']!=metadata[row['state']]['exon']:raise ValueError('Source identity changed')
        if type(row['source_row']) is not int or row['source_row'] in source_index:raise ValueError('Duplicate source row')
        source_index[row['source_row']]=key;by[key].append(row['source_row'])
    if len(source_rows)!=396 or any(len(by[key,cell])!=3 for key in metadata for cell in CELLS):raise ValueError('Metadata multiplicity')
    for row in states:
        key=(row['state'],row['cell'])
        if row['state'] not in metadata or row['cell'] not in CELLS or key in interval_by:raise ValueError('State-cell identity or duplication')
        if sorted(row['source_rows'])!=sorted(by[key]) or type(row['usable_rows']) is not int or not 0<=row['usable_rows']<=3 or type(row['complete']) is not bool or row['complete']!=(row['usable_rows']==3):raise ValueError('Interval metadata semantics')
        if not isinstance(row['lower'],str) or (row['upper'] is not None and not isinstance(row['upper'],str)):raise ValueError('Exact interval strings required')
        lo=F(row['lower']);up=None if row['upper'] is None else F(row['upper'])
        if row['complete']:
            if lo<=0 or up is None or up<lo:raise ValueError('Complete positive interval invalid')
        elif lo!=0 or up is not None:raise ValueError('Unknown must remain (0,infinity)')
        interval_by[key]=(lo,up)
    if len(interval_by)!=132:raise ValueError('Interval coverage changed')
    previous={}
    if [F(s['gamma']) for s in sensitivities]!=GAMMAS:raise ValueError('Sensitivity grid changed')
    for s in sensitivities:
        if [x['exon'] for x in s['exon_shared_order_results']]!=EXONS:raise ValueError('Existing order coverage')
        for x in s['exon_shared_order_results']:
            if x['status'] not in ['feasible','infeasible']:raise ValueError('Existing order unresolved')
            previous[str(F(s['gamma'])),x['exon']]=x['status']
    records=[{'cell':r['cell'],'state':r['id'],'exon':r['exon']} for r in expected]
    intervals=[interval_by[r['state'],r['cell']] for r in records]
    for X in matrices['matrices'].values():validate_design(X,records,intervals)
    return records,intervals,matrices['matrices'],previous

def analyze(source,matrix_path,mapping_path,checkpoint):
    records,base,matrices,previous=load_contract(source,matrix_path,mapping_path)
    results=[];statusmap={};unknown=sum(u is None for l,u in base)
    for gamma in GAMMAS:
        intervals=[(l/gamma,None if u is None else u*gamma) for l,u in base]
        for exon in EXONS:
            for model in MODELS:
                arms=['both'] if model==MODELS[0] else ['both']+CELLS
                for arm in arms:
                    inds=[i for i,r in enumerate(records) if r['exon']==exon and (arm=='both' or r['cell']==arm)]
                    X=[matrices[model][i] for i in inds];rs=[records[i] for i in inds];ints=[intervals[i] for i in inds]
                    item={'gamma':str(gamma),'exon':exon,'model':model,'arm':arm,'row_indices':inds,
                      'unknown_intervals':sum(u is None for l,u in ints),'state_cell_count':len(inds)}
                    try:item.update(solve(X,rs,ints))
                    except ArithmeticError as exc:
                        item.update(status='unresolved',error=str(exc));results.append(item);checkpoint(results)
                        raise ArithmeticError('Unresolved certificate: '+str(item)) from exc
                    results.append(item);checkpoint(results)
                    if model==MODELS[2] and item['status']!='feasible':raise RuntimeError('Saturated pair model failed')
                    if model==MODELS[1] and exon!=17 and item['status']!='feasible':raise RuntimeError('Full-rank additive model failed')
                    if model==MODELS[0] and item['status']=='feasible' and previous[str(gamma),exon]=='infeasible':raise RuntimeError('Additive model contradicts weaker shared order')
                    statusmap[str(gamma),exon,model,arm]=item['status']
    if len(results)!=140 or sum(r['arm']=='both' for r in results)!=60:raise RuntimeError('Declared coverage failed')
    for exon in EXONS:
        for model in MODELS:
            for arm in (['both'] if model==MODELS[0] else ['both']+CELLS):
                series=[statusmap[str(g),exon,model,arm] for g in GAMMAS]
                if any(series[i]=='feasible' and series[j]=='infeasible' for i in range(4) for j in range(i+1,4)):raise RuntimeError('Gamma nesting failure')
        for g in GAMMAS:
            joint=[statusmap[str(g),exon,m,'both'] for m in MODELS]
            if any(joint[i]=='feasible' and joint[j]=='infeasible' for i in range(3) for j in range(i+1,3)):raise RuntimeError('Model nesting failure')
            for model in MODELS[1:]:
                if (statusmap[str(g),exon,model,'both']=='feasible')!=all(statusmap[str(g),exon,model,c]=='feasible' for c in CELLS):raise RuntimeError('Context block conjunction failed')
    summary=[]
    for g in GAMMAS:
        summary.append({'gamma':str(g),'joint_infeasible_exons':{m:[e for e in EXONS if statusmap[str(g),e,m,'both']=='infeasible'] for m in MODELS},
          'context_additive_infeasible':{c:[e for e in EXONS if statusmap[str(g),e,MODELS[1],c]=='infeasible'] for c in CELLS},
          'existing_unrestricted_common_order_infeasible_exons':[e for e in EXONS if previous[str(g),e]=='infeasible']})
    return {'status':'COMPLETE_EXPLORATORY_CERTIFIED_AUDIT','state_cells':132,'unknown_state_cells':unknown,
      'joint_cases':60,'individual_cases':80,'gamma_and_model_nesting_passed':True,'summaries':summary,'results':results,
      'scope':'Exploratory reuse of exposed local-state proxy intervals. Recorded descriptors, not complete haplotypes. No individual interaction identification, calibrated confidence, causal attribution or new validation.'}

def runtime():
    import numpy,scipy
    return {'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__}

def main():
    if not __debug__:raise RuntimeError('Optimized execution is prohibited before input access')
    ap=argparse.ArgumentParser();ap.add_argument('--freeze',type=Path,required=True);ap.add_argument('--authorization',type=Path,required=True)
    ap.add_argument('--b1',type=Path,required=True);ap.add_argument('--matrices',type=Path,required=True);ap.add_argument('--mapping',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.out.exists():raise FileExistsError(args.out)
    auth=json.loads(args.authorization.read_text());freeze=json.loads(args.freeze.read_text())
    if auth.get('permit_exploratory_fit') is not True or auth.get('freeze_sha256')!=sha(args.freeze):raise ValueError('Matching exploratory authorization required')
    here=Path(__file__).resolve().parent
    for name,h in freeze['code_sha256'].items():
        if sha(here/name)!=h:raise ValueError('Frozen implementation changed: '+name)
    paths={'b1':args.b1,'matrices':args.matrices,'mapping':args.mapping}
    for role,p in paths.items():
        if sha(p)!=freeze['input_sha256'][role]:raise ValueError('Input hash mismatch: '+role)
    if runtime()!=freeze['runtime']:raise RuntimeError('Runtime differs from reviewed freeze')
    args.out.mkdir(parents=True);started=time.perf_counter()
    save(args.out/'invocation.json',{'freeze_sha256':sha(args.freeze),'authorization_sha256':sha(args.authorization),'input_sha256':freeze['input_sha256'],'runtime':runtime(),'command_arguments':{k:str(v) for k,v in vars(args).items()}})
    def checkpoint(rs):save(args.out/'progress.json',{'status':'IN_PROGRESS','completed_cases':len(rs),'results':rs})
    try:
        result=analyze(args.b1,args.matrices,args.mapping,checkpoint)
        for role,p in paths.items():
            if sha(p)!=freeze['input_sha256'][role]:raise RuntimeError('Input changed during execution')
        save(args.out/'result.json',result)
        save(args.out/'receipt.json',{'status':result['status'],'result_sha256':sha(args.out/'result.json'),'freeze_sha256':sha(args.freeze),'input_sha256':freeze['input_sha256'],'runtime':runtime(),'elapsed_seconds':time.perf_counter()-started,'scope':result['scope']})
        print(json.dumps({'status':result['status'],'summaries':result['summaries']}))
    except Exception as exc:
        save(args.out/'failure.json',{'status':'STOPPED_UNRESOLVED','type':type(exc).__name__,'error':str(exc),'freeze_sha256':sha(args.freeze),'elapsed_seconds':time.perf_counter()-started});raise

if __name__=='__main__':main()
