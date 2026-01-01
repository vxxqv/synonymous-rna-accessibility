"""Independent small-case references for the held synthetic-only model kernel."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import importlib.util,hashlib,json,random

if not __debug__:raise RuntimeError('Optimized Python prohibited')
p=Path('work/research_cycle8/comt_independent_review')
dest=p/'independent_kernel_checks.json'
if dest.exists():raise FileExistsError(dest)
src=Path('work/research_cycle8/comt_pilot/proxy_models.py')
assert hashlib.sha256(src.read_bytes()).hexdigest()=='8811af9750e813c798122c7e936b6279a0f3b0e53fea1924cb2d633162a34093'
spec=importlib.util.spec_from_file_location('reviewed_kernel',src);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
checks=0

def exact_iso_reference(x,y,w):
    # Enumerate every contiguous partition after exact ties are consolidated.
    # Choose nondecreasing block means with minimum weighted SSE.
    xs=sorted(set(x));yy=[];ww=[]
    for a in xs:
        idx=[i for i,z in enumerate(x) if z==a]
        mass=sum(w[i] for i in idx);ww.append(F(mass));yy.append(sum(F(w[i])*y[i] for i in idx)/mass)
    best=None;solution=None
    for cuts in product([False,True],repeat=len(xs)-1):
        bounds=[0]+[i+1 for i,c in enumerate(cuts) if c]+[len(xs)]
        vals=[]
        for a,b in zip(bounds,bounds[1:]):
            v=sum(ww[i]*yy[i] for i in range(a,b))/sum(ww[a:b]);vals.extend([v]*(b-a))
        if any(a>b for a,b in zip(vals,vals[1:])):continue
        sse=sum(ww[i]*(yy[i]-vals[i])**2 for i in range(len(xs)))
        if best is None or sse<best:best,solution=sse,vals
    return xs,solution

rng=random.Random(731908)
for k in range(80):
    x=[rng.randrange(-2,3) for _ in range(6)]
    z=[rng.randrange(-2,3) for _ in range(6)]
    y=[rng.randrange(-4,5) for _ in range(6)]
    w=[rng.randrange(1,4) for _ in range(6)]
    xs,expected=exact_iso_reference(x,y,w);got=m.fit_pava(x,y,w)
    assert got['x']==xs and got['y']==expected;checks+=1
    # KKT sufficiency for the convex constrained weighted LS problem.
    fit=m.fit_linear(x,y,w,z);a,b,c=fit['rna_slope'],fit['load_slope'],fit['intercept']
    residual=[F(yi)-c-a*xi-b*zi for xi,yi,zi in zip(x,y,z)]
    assert a>=0 and sum(t*r for t,r in zip(w,residual))==0
    assert sum(t*zi*r for t,zi,r in zip(w,z,residual))==0
    cross=sum(t*xi*r for t,xi,r in zip(w,x,residual))
    assert (a>0 and cross==0) or (a==0 and cross<=0);checks+=1

# Model and score checks with unequal numbers of alleles per site.
weights=m.site_weights(['s1','s1','s1','s2'])
assert weights==[F(1,3),F(1,3),F(1,3),F(1)]
assert m.site_mae([0]*4,[1,1,1,3],['s1','s1','s1','s2'])==2;checks+=2
roi=m.fit_roi_medians([0,2,100,102],[1]*4,['ROI_1']*2+['ROI_2']*2)
assert roi=={'ROI_1':1,'ROI_2':101};checks+=1
assert m.expansion_criterion(F(9),[10]*4,100,40)['pass']
assert not m.expansion_criterion(F(9001,1000),[10]*4,100,40)['pass']
assert not m.expansion_criterion(0,[0,10,10,10],125,53)['pass'];checks+=3

runtime=Path('work/research_cycle8/comt_pilot/runtime')
manifest=json.loads((runtime/'runtime_manifest.json').read_text())
filechecks=[]
for relative,wanted in manifest['package_files'].items():
    f=runtime/'library'/relative
    # Manifest paths are data, never shell commands.
    assert f.resolve().is_relative_to((runtime/'library').resolve())
    got=hashlib.sha256(f.read_bytes()).hexdigest()
    if got!=wanted:filechecks.append(relative)
assert not filechecks
assert hashlib.sha256(Path(manifest['Rscript']).read_bytes()).hexdigest()==manifest['Rscript_sha256']
out={'status':'PASS_BOUNDED_INDEPENDENT_REFERENCES','kernel_checks':checks,
     'PAVA_reference':'Exhaustive contiguous-partition minimum weighted SSE with exact tie consolidation',
     'linear_reference':'Exact KKT sufficient conditions of constrained convex least squares',
     'runtime_package_file_hashes_checked':len(manifest['package_files']),
     'runtime_hash_mismatches':filechecks,'Rscript_hash_matches':True,
     'runtime_smoke_not_reexecuted':True,'raw_reader_or_staged_driver_verified':False,
     'no_source_count_or_phenotype_values_read':True,
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':checks,'runtime_files':len(manifest['package_files']),'receipt_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}))
