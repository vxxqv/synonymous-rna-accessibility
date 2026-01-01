"""Exact minimal positive dependencies of size <=3; exposed FAS application."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import argparse, hashlib, json, math

PIN='3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466'

def positive_weights(rows):
    n=len(rows)
    if n==1:
        return [Q(1)] if not any(rows[0]) else None
    if n==2:
        a,b=rows
        k=next((i for i,v in enumerate(a) if v),None)
        if k is None:return None # nonminimal zero vector case
        r=-Q(b[k],a[k])
        if r>0 and all(r*x+y==0 for x,y in zip(a,b)):
            return [r/(r+1),1/(r+1)]
        return None
    if n!=3:raise ValueError('Scope is one through three rows')
    # A minimal three-row dependence has rank2. Any nonzero 2x2 minor
    # provides its unique null direction; all remaining columns are checked.
    for i,j in combinations(range(len(rows[0])),2):
        a=[r[i] for r in rows];b=[r[j] for r in rows]
        w=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
        if not any(w):continue
        if all(x<0 for x in w):w=[-x for x in w]
        if all(x>0 for x in w) and all(sum(w[k]*rows[k][c] for k in range(3))==0 for c in range(len(rows[0]))):
            return [Q(x,sum(w)) for x in w]
        return None
    return None

def fixtures():
    assert positive_weights([[0,0]])==[Q(1)]
    assert positive_weights([[1,0]]) is None
    assert positive_weights([[1,2],[-2,-4]])==[Q(2,3),Q(1,3)]
    assert positive_weights([[1,2],[2,4]]) is None
    assert positive_weights([[1,0],[0,1],[-1,-1]])==[Q(1,3)]*3
    assert positive_weights([[1,0],[0,1],[1,1]]) is None
    assert positive_weights([[1,0,0],[0,1,0],[-1,-1,1]]) is None
    assert positive_weights([[1,0],[2,0],[-1,0]]) is None # nonminimal triple
    assert not (Q(2)*Q(1) < Q(2)) # strict comparison lost at equality
    return 9

def analyze(d):
    v=d['vertices_data']
    if len(v)!=72 or [r['vertex'] for r in v]!=list(range(72)):raise ValueError('Domain changed')
    s=[r['codon_states'] for r in v]
    levels=[sorted({r[k] for r in s}) for k in range(len(s[0]))]
    x=[[1]+[int(r[k]==a) for k,ls in enumerate(levels) for a in ls[1:]] for r in s]
    if len(x[0])!=8:raise ValueError('Design changed')
    l=[Q(r['lower']) for r in v];u=[Q(r['upper']) for r in v]
    if l[46]!=1 or u[46]!=1 or any(a<=0 or a>b for a,b in zip(l,u)):raise ValueError('Interval domain changed')
    edges=[]
    for i in range(72):
        for j in range(72):
            if u[i]<l[j]:
                a=int(i!=46)+int(j!=46)
                edges.append({'i':i,'j':j,'row':[b-a for a,b in zip(x[i],x[j])], 'loss_t':str((l[j]/u[i])**(2//a))})
    found=[];smaller=[];checks={}
    for n in (1,2,3):
        total=0;skipped=0
        for ids in combinations(range(len(edges)),n):
            total+=1
            if any(set(z)<=set(ids) for z in smaller):
                skipped+=1;continue
            w=positive_weights([edges[k]['row'] for k in ids])
            if w:
                if any(z<=0 for z in w) or sum(w)!=1:raise ArithmeticError('Invalid positive weights')
                assert all(sum(w[q]*edges[k]['row'][c] for q,k in enumerate(ids))==0 for c in range(8))
                t=min(Q(edges[k]['loss_t']) for k in ids)
                found.append({'edge_indices':list(ids),'support_size':n,'edges':[edges[k] for k in ids], 'weights':list(map(str,w)), 'vertices':sorted({a for k in ids for a in (edges[k]['i'],edges[k]['j'])}), 'persistence_t_exclusive':str(t)})
                smaller.append(ids)
        checks[str(n)]={'combinations':total,'nonminimal_skipped':skipped,'certificates':sum(r['support_size']==n for r in found)}
    minimum=min((r['support_size'] for r in found),default=None)
    best=max((Q(r['persistence_t_exclusive']) for r in found),default=None)
    common=set(found[0]['vertices']) if found else set()
    for r in found:common.intersection_update(r['vertices'])
    return {'scope':'Exposed exploratory minimum-comparison audit; all minimal supports of size at most3, not all minimal supports.', 'states':72,'design_columns':8,'forced_orders':len(edges),'enumeration':checks,'minimum_comparisons':minimum,'certificates':found,'all_forced_orders':edges,'common_vertices_across_small_certificates':sorted(common),'strongest_small_certificate_t':None if best is None else str(best),'global_boundary_previously_verified':'9128860950384/8609574801263', 'small_certificate_reaches_global_boundary':best==Q(9128860950384,8609574801263)}

def main():
    if not __debug__:raise RuntimeError('No optimized mode')
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('No overwrite')
    n=fixtures();raw=a.input.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PIN:raise ValueError('Changed source')
    r=analyze(json.loads(raw));r['synthetic_assertions']=n;r['source_sha256']=PIN
    a.out.mkdir(parents=True);(a.out/'result.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    m={'analysis_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256(Path(__file__).with_name('protocol.md').read_bytes()).hexdigest(),'result_sha256':hashlib.sha256((a.out/'result.json').read_bytes()).hexdigest()}
    (a.out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:r[k] for k in ('forced_orders','enumeration','minimum_comparisons','common_vertices_across_small_certificates','strongest_small_certificate_t','small_certificate_reaches_global_boundary')}))

if __name__=='__main__':main()
