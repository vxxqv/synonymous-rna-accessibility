"""Synthetic-only reference and guard checks; no observed B1 result reads."""
import copy, hashlib, itertools, json, random, subprocess, sys, tempfile
from fractions import Fraction as F
from pathlib import Path
import runner as r

HERE=Path(__file__).resolve().parent
N=0
def require(value,message):
    global N
    if not value:raise AssertionError(message)
    N+=1

def records(n,cell='K562'):
    return [{'state':str(i),'cell':cell,'exon':17} for i in range(n)]

def expect_error(fn,typ=Exception):
    try:fn()
    except typ:return True
    raise AssertionError('Expected refusal')

def analytic_rectangle(intervals):
    """Enumerate 24 strict latent orders using an independent gap-cone test.

    For s=(1,-1,-1,1), sum s_i z_i=0 becomes a linear form in
    three strictly positive adjacent gaps. Both coefficient signs are necessary
    and sufficient. A feasible strict-inequality system can be generically
    perturbed within the additive subspace to avoid latent ties.
    """
    signs=[1,-1,-1,1]
    for p in itertools.permutations(range(4)):
        pos={v:i for i,v in enumerate(p)}
        if any(up is not None and up<intervals[j][0] and pos[i]>=pos[j] for i,(lo,up) in enumerate(intervals) for j in range(4)):continue
        coeff=[sum(signs[p[j]] for j in range(k+1,4)) for k in range(3)]
        if min(coeff)<0<max(coeff):return 'feasible'
    return 'infeasible'

def main():
    X=[[1,0,0],[1,0,1],[1,1,0],[1,1,1]];rs=records(4)
    interval_types=[(F(1),F(1)),(F(1),F(2)),(F(2),F(2)),(F(0),None)]
    for intervals in itertools.product(interval_types,repeat=4):
        ans=r.solve(X,rs,intervals)
        require(ans['status']==analytic_rectangle(intervals),'Independent finite gap-cone disagreement')
        require(r.verify_certificate(X,rs,intervals,ans),'Exact certificate')
    for p in itertools.permutations(range(4)):
        levels={v:F(i+1) for i,v in enumerate(p)};ints=[(levels[i],levels[i]) for i in range(4)]
        ans=r.solve(X,rs,ints);require(ans['status']==analytic_rectangle(ints),'Strict-order oracle')
    rng=random.Random(3108);checker=[(F(2),F(2)),(F(1),F(1)),(F(1),F(1)),(F(2),F(2))]
    for k in range(12):
        ix=list(range(4));jx=list(range(3));rng.shuffle(ix);rng.shuffle(jx)
        xx=[[X[i][j] for j in jx] for i in ix];rr=[rs[i] for i in ix];ii=[checker[i] for i in ix]
        ans=r.solve(xx,rr,ii);require(ans['status']=='infeasible','Permutation changed result')
    # Column complements preserve span, even with arbitrary affine recoding.
    rec=[[row[0],row[0]-row[1],row[0]-row[2]] for row in X]
    require(r.solve(rec,rs,checker)['status']=='infeasible','Coding changed status')
    sat=[[int(i==j) for j in range(4)] for i in range(4)]
    good=r.solve(sat,rs,checker);require(good['status']=='feasible','Saturated field rejected')
    bad=copy.deepcopy(good);bad['responses'][0]='0';require(expect_error(lambda:r.verify_certificate(sat,rs,checker,bad),ArithmeticError),'Bad completion accepted')
    dual=r.solve(X,rs,checker);bad=copy.deepcopy(dual);bad['dual'][0]['weight']='0';require(expect_error(lambda:r.verify_certificate(X,rs,checker,bad),ArithmeticError),'Bad dual accepted')
    same=[[1],[1],[1]];sameints=[(F(1),F(2)),(F(2),F(3)),(F(0),None)]
    tied=r.solve(same,records(3),sameints);require(tied['responses']==['2','2','2'],'Tied unknown completion')
    unknown=r.solve([[1],[1]],records(2),[(F(0),None)]*2);require(unknown['responses']==['1','1'],'Unknown positive domain')
    require(expect_error(lambda:r.solve([[True]],records(1),[(F(1),F(1))]),ValueError),'Boolean matrix')
    require(expect_error(lambda:r.solve([[1]],records(1),[(F(0),F(0))]),ValueError),'Empty positive interval')
    # Independent cell responses: cross-cell order is never accidentally imposed.
    shared=[[1,0],[1,1],[1,0],[1,1]];ctx=records(2,'K562')+records(2,'KBM7')
    ints=[(F(1),F(1)),(F(2),F(2)),(F(2),F(2)),(F(1),F(1))]
    require(r.solve(shared,ctx,ints)['status']=='infeasible','Shared reversal should fail')
    separate=[[1,0,0,0],[1,1,0,0],[0,0,1,0],[0,0,1,1]]
    require(r.solve(separate,ctx,ints)['status']=='feasible','Separate contexts rejected')
    # Full metadata universe, invented constant intervals, forbidden source fields.
    mp=HERE.parent/'inventory_v1/nested_matrices.json';mapfile=HERE.parent/'exact_v1/mapping.json'
    mapping=json.loads(mapfile.read_text());meta={x['id']:x for x in mapping};states=[];source=[];j=100
    for key in sorted(meta):
        for cell in r.CELLS:
            nums=[]
            for z in range(3):
                nums.append(j);source.append({'source_row':j,'state':key,'cell':cell,'exon':meta[key]['exon'],
                  'DNA':'POISON_UNSELECTED','RNA':'POISON_UNSELECTED','ratio':'POISON_UNSELECTED','source_score':'POISON_UNSELECTED'});j+=1
            states.append({'state':key,'cell':cell,'usable_rows':3,'complete':True,'lower':'1','upper':'2','source_rows':nums,'median':'POISON_UNSELECTED'})
    sens=[{'gamma':str(g),'group_cell_results':'POISON_UNSELECTED','exon_shared_order_results':[{'exon':e,'status':'feasible','edges':'POISON_UNSELECTED'} for e in r.EXONS]} for g in r.GAMMAS]
    toy={'status':'COMPLETE_DESCRIPTIVE_AUDIT','source_rows':source,'state_cells':states,'sensitivities':sens}
    with tempfile.TemporaryDirectory(prefix='synthetic_',dir=HERE) as td:
        d=Path(td);p=d/'toy.json';r.save(p,toy)
        original=r.json.loads
        def forbid(s,*a,**kw):
            if 'POISON_UNSELECTED' in s:raise AssertionError('Unselected field decoded')
            return original(s,*a,**kw)
        r.json.loads=forbid
        try:read=r.parse_b1(p)
        finally:r.json.loads=original
        require(len(read[0])==132 and len(read[1])==396,'Selective coverage')
        allout=r.analyze(p,mp,mapfile,lambda _:None)
        require(len(allout['results'])==140 and all(x['status']=='feasible' for x in allout['results']),'Full synthetic protocol coverage')
        broken=copy.deepcopy(toy);broken['state_cells'].append(broken['state_cells'][0]);bp=d/'duplicate.json';r.save(bp,broken)
        require(expect_error(lambda:r.load_contract(bp,mp,mapfile),ValueError),'Duplicate state accepted')
        broken=copy.deepcopy(toy);broken['state_cells'][0]['lower']=1;r.save(bp,broken)
        require(expect_error(lambda:r.load_contract(bp,mp,mapfile),ValueError),'Inexact bound type accepted')
        # Invalid hash must fail before decoding an invented poisoned source.
        fp=d/'freeze.json';auth=d/'auth.json';poison=d/'poison.json';poison.write_text('INVALID_SOURCE_MUST_NOT_BE_DECODED')
        freeze={'code_sha256':{'runner.py':r.sha(HERE/'runner.py')},'input_sha256':{'b1':'0'*64,'matrices':r.sha(mp),'mapping':r.sha(mapfile)},'runtime':r.runtime()}
        r.save(fp,freeze);r.save(auth,{'permit_exploratory_fit':True,'freeze_sha256':r.sha(fp)})
        out=d/'refused'
        args=[sys.executable,str(HERE/'runner.py'),'--freeze',str(fp),'--authorization',str(auth),'--b1',str(poison),'--matrices',str(mp),'--mapping',str(mapfile),'--out',str(out)]
        proc=subprocess.run(args,capture_output=True,text=True)
        require(proc.returncode!=0 and 'Input hash mismatch: b1' in proc.stderr and not out.exists(),'Hash guard order')
        for opt in ['-O','-OO']:
            proc=subprocess.run([sys.executable,opt,str(HERE/'runner.py')],capture_output=True,text=True)
            require(proc.returncode!=0 and 'Optimized execution is prohibited' in proc.stderr,'Optimization guard')
        out.mkdir();proc=subprocess.run(args,capture_output=True,text=True)
        require(proc.returncode!=0 and 'FileExistsError' in proc.stderr,'Overwrite guard')
    receipt={'status':'PASS','checks':N,'independent_rectangle_interval_cases':256,'strict_order_cases':24,
      'full_synthetic_protocol_records':140,'outcome_values_read':False,
      'runner_sha256':r.sha(HERE/'runner.py'),'test_source_sha256':r.sha(Path(__file__)),'runtime':r.runtime(),
      'limits':'Independent finite gap-cone oracle within the same authored implementation lane. Root will separately verify observed certificates.'}
    if (HERE/'synthetic_receipt.json').exists():raise FileExistsError('Preserve previous test receipt')
    r.save(HERE/'synthetic_receipt.json',receipt);print(json.dumps(receipt))

if __name__=='__main__':main()
