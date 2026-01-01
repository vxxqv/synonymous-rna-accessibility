"""Exploratory exact robustness of a fixed ordinal model, with a fixed reference."""
from fractions import Fraction as F
from pathlib import Path
import argparse, gzip, hashlib, importlib.util, json, math, sys

def load_core(path):
    spec=importlib.util.spec_from_file_location('frozen_ordinal',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def fields(lower,upper,mask,value,metric):
    v=F(value)
    if metric=='additive':
        return ([max(F(0),l-a*v) for l,a in zip(lower,mask)],
                [None if u is None else u+a*v for u,a in zip(upper,mask)])
    if metric!='multiplicative_squared' or v<1:raise ValueError('Invalid metric/parameter')
    # Squaring is strictly increasing on nonnegative scores.
    return ([l*l/(v**a) for l,a in zip(lower,mask)],
            [None if u is None else u*u*(v**a) for u,a in zip(upper,mask)])

def breakpoints(lower,upper,mask,metric):
    values={F(0) if metric=='additive' else F(1)}
    thresholds={}
    for i,u in enumerate(upper):
        if u is None:continue
        for j,l in enumerate(lower):
            if u>=l:continue
            a=mask[i]+mask[j]
            if a==0 or (metric=='multiplicative_squared' and u==0):
                thresholds[(i,j)]=None
            else:
                t=(l-u)/a if metric=='additive' else (l/u)**(2//a)
                thresholds[(i,j)]=t;values.add(t)
    return sorted(values),thresholds

def sqrt_enclosure(value,places=12):
    v=F(value);scale=10**places
    k=math.isqrt((v.numerator*scale*scale)//v.denominator)
    a,b=F(k,scale),F(k+1,scale)
    if a*a==v:b=a
    if not a*a<=v<=b*b:raise ArithmeticError('Root enclosure failed')
    def decimal(integer):
        return str(integer//scale)+'.'+str(integer%scale).zfill(places)
    return {'lower':str(a),'upper':str(b),'decimal_lower_exact':decimal(k),
            'decimal_upper_exact':decimal(k if a==b else k+1)}

def radius(core,x,lower,upper,mask,metric):
    lower=list(map(F,lower));upper=[None if u is None else F(u) for u in upper]
    if not lower or len(lower)!=len(upper) or len(mask)!=len(lower):raise ValueError('Shape mismatch')
    if any(l<0 or (u is not None and u<l) for l,u in zip(lower,upper)):raise ValueError('Nonnegative nonempty intervals required')
    if any(type(a) is not int or a not in (0,1) for a in mask):raise ValueError('Integer binary perturbation mask required')
    candidates,thresholds=breakpoints(lower,upper,mask,metric)
    tested={}
    def test(k):
        if k not in tested:
            lo,hi=fields(lower,upper,mask,candidates[k],metric)
            tested[k]=core.order_certificate(x,lo,hi)
        return tested[k]
    if test(0)['status']=='feasible':at=0;before=None
    elif test(len(candidates)-1)['status']=='infeasible':
        cert=tested[len(candidates)-1]
        if any(thresholds[(d['lower_vertex'],d['higher_vertex'])] is not None for d in cert['dual']):
            raise ArithmeticError('Infinite radius needs permanent contradiction')
        return {'metric':metric,'status':'infinite','permanent_dual':cert,
                'breakpoints':list(map(str,candidates)),'tested':{str(k):v for k,v in sorted(tested.items())}}
    else:
        a,b=0,len(candidates)-1
        while b-a>1:
            mid=(a+b)//2
            if test(mid)['status']=='feasible':b=mid
            else:a=mid
        at=b;before=a
    value=candidates[at];feasible=test(at)
    lower_certificate=None
    if before is not None:
        lower_certificate=test(before)
        for d in lower_certificate['dual']:
            t=thresholds[(d['lower_vertex'],d['higher_vertex'])]
            if t is not None and t<value:raise ArithmeticError('Dual does not persist to sharp boundary')
    lo,hi=fields(lower,upper,mask,value,metric)
    result={'metric':metric,'status':'finite','parameter_exact':str(value),
        'parameter_meaning':'epsilon' if metric=='additive' else 'gamma_squared',
        'sharp_lower_certificate':lower_certificate,'feasible_at_boundary':feasible,
        'boundary_lower':list(map(str,lo)),'boundary_upper':[None if u is None else str(u) for u in hi],
        'response_meaning':'original units' if metric=='additive' else 'square roots of returned rational responses',
        'breakpoints':list(map(str,candidates)), 'tested':{str(k):v for k,v in sorted(tested.items())}}
    if metric!='additive':result['gamma_enclosure']=sqrt_enclosure(value)
    return result

def synthetic(core):
    x=[[1,0,0],[1,1,0],[1,0,1],[1,1,1]];mask=[0,1,1,1]
    point=list(map(F,['1','1/4','1/4','1']))
    a=radius(core,x,point,point,mask,'additive')
    m=radius(core,x,point,point,mask,'multiplicative_squared')
    if a['parameter_exact']!='3/8' or m['parameter_exact']!='4':raise ArithmeticError('Reciprocal fixture')
    z=list(map(F,[1,0,0,1]))
    if radius(core,x,z,z,mask,'multiplicative_squared')['status']!='infinite':raise ArithmeticError('Zero upper permanent fixture')
    const=[F(1)]*4
    if radius(core,x,const,const,mask,'additive')['parameter_exact']!='0':raise ArithmeticError('Constant feasible fixture')
    if radius(core,x,const,const,mask,'multiplicative_squared')['parameter_exact']!='1':raise ArithmeticError('Constant multiplicative fixture')
    lower=[F(0),F(1),F(0),F(1)];upper=[None,F(1),F(0),F(1)]
    unknown=radius(core,x,lower,upper,mask,'additive')
    if unknown['parameter_exact']!='0':raise ArithmeticError('Unknown feasible fixture')
    fixed=radius(core,[[1],[1]],[F(0),F(1)],[F(0),F(1)],[0,0],'additive')
    if fixed['status']!='infinite':raise ArithmeticError('Fixed contradictory observations')
    # Anchored one-factor duplicate forces a linear gamma breakpoint.
    dup=radius(core,[[1],[1]],[F(1),F(1,2)],[F(1),F(1,2)],[0,1],'multiplicative_squared')
    if dup['parameter_exact']!='4':raise ArithmeticError('Anchor exponent fixture')
    roots=[sqrt_enclosure(q) for q in [F(1),F(2),F(4),F(1000001,1000000)]]
    return {'status':'passed','fixtures':8,'root_enclosures':roots,'reciprocal_additive':a,'reciprocal_multiplicative':m}

def analyze(core,d):
    vertices=d['vertices_data']
    if len(vertices)!=72 or [v['vertex'] for v in vertices]!=list(range(72)):
        raise ValueError('Requires the exact 72-vertex domain')
    if any(v['upper'] is None or len(v['scenario_ratios'])!=27 or any(s is None for s in v['scenario_ratios']) for v in vertices):
        raise ValueError('Observed adapter requires complete finite input; it must not drop unknown states')
    x=core.features([v['codon_states'] for v in vertices])
    lo=[F(v['lower']) for v in vertices];hi=[None if v['upper'] is None else F(v['upper']) for v in vertices]
    if len(vertices)!=72 or vertices[46]['lower']!='1' or vertices[46]['upper']!='1':raise ValueError('Frozen FAS domain/reference mismatch')
    mask=[0 if i==46 else 1 for i in range(72)]
    metrics=['additive','multiplicative_squared']
    full={m:radius(core,x,lo,hi,mask,m) for m in metrics}
    deletions=[]
    for omit in range(72):
        ids=[i for i in range(72) if i!=omit]
        cert=core.order_certificate([x[i] for i in ids],[lo[i] for i in ids],[hi[i] for i in ids])
        deletions.append({'omitted_vertex':omit,'retained_vertex_ids':ids,'certificate':cert,
                          'note':'Reference still used in original normalization even when its observation row is omitted.'})
    prep=[]
    for kind,n in [('output',9),('input',3)]:
        for omitted in range(n):
            cols=[k for k in range(27) if (k//3 if kind=='output' else k%3)!=omitted]
            vals=[[F(v['scenario_ratios'][k]) for k in cols] for v in vertices]
            l=[min(v) for v in vals];u=[max(v) for v in vals]
            if any(a>b or c<e for a,b,c,e in zip(lo,l,hi,u)):raise ArithmeticError('Deletion envelope must nest')
            res={m:radius(core,x,l,u,mask,m) for m in metrics}
            for m in metrics:
                if res[m]['status']=='finite' and F(res[m]['parameter_exact'])<F(full[m]['parameter_exact']):raise ArithmeticError('Nested-set radius monotonicity')
            prep.append({'kind':kind,'omitted_zero_based':omitted,'retained_scenarios':cols,'results':res})
    return {'scope':'Exploratory post-primary exact sensitivity; not fresh validation or biological calibration.',
        'anchor_vertex':46,'full':full,'all_single_vertex_deletions':deletions,
        'all_single_preparation_deletions':prep,
        'preparation_deletion_interpretation':'Removing scenarios narrows the envelope and cannot repair an exclusion; these are influence diagnostics, not independent robustness evidence.'}

def main():
    if not __debug__:raise RuntimeError('Optimized Python is prohibited')
    p=argparse.ArgumentParser();p.add_argument('--core',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--result',type=Path);p.add_argument('--freeze',type=Path);args=p.parse_args()
    if args.out.exists():raise FileExistsError('Refusing an existing output')
    if args.result:
        if not args.freeze:raise ValueError('Observed successor requires a frozen contract')
        frozen=json.loads(args.freeze.read_text())
        for name,sha in frozen['files'].items():
            if hashlib.sha256((args.freeze.parent/name).read_bytes()).hexdigest()!=sha:raise ValueError('Changed frozen follow-up file: '+name)
        raw=args.result.read_bytes()
        if args.result.suffix=='.gz':raw=gzip.decompress(raw)
        if hashlib.sha256(raw).hexdigest()!=frozen['primary_result_sha256']:raise ValueError('Changed primary result')
        if hashlib.sha256(args.core.read_bytes()).hexdigest()!=frozen['ordinal_core_sha256']:raise ValueError('Changed ordinal core')
        core=load_core(args.core)
        d=analyze(core,json.loads(raw));d['source_sha256']=hashlib.sha256(raw).hexdigest()
    else:d=synthetic(load_core(args.core))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'completed','output':str(args.out)}))

if __name__=='__main__':main()
