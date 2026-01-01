"""Independent exact check of the new boundary witnesses; no solver import."""
import json, hashlib
from pathlib import Path
from fractions import Fraction as Q
from sympy import Matrix

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if not __debug__: raise RuntimeError('Assertions required')
    base=Path(__file__).resolve().parents[2]
    files={'result':base/'context_transport_sharp_boundary/actual_v1/result.json',
           'mapping':base/'context_design_audit/exact_v1/mapping.json',
           'source':base/'cprime_results/primary_v1/result.json'}
    pins={'result':'3477a640e9b16220498edbdf92ee7a2fd46292dc715c1fb58136344e93293afc',
          'mapping':'a2f428f5308cb78a827bfbd23f7d7756b98980e3585f0ccf6a733859c2388436',
          'source':'eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'}
    assert all(sha(files[k])==v for k,v in pins.items())
    r,mp,src=[json.loads(files[k].read_text()) for k in ('result','mapping','source')]
    mp.sort(key=lambda z:z['id']); es=sorted({z['exon'] for z in mp})
    aa=sorted({a for z in mp for a in z['alleles']})
    x=[[int(z['exon']==e) for e in es]+[int(a in z['alleles']) for a in aa] for z in mp]
    rows=[dict(z,cell=c) for c in ('K562','KBM7') for z in mp]
    lu={(z['state'],z['cell']):(Q(z['lower']), Q(z['upper'])) for z in src['state_cells']}
    totals={'primal_cases':0,'forced_comparisons':0,'responses':0,'isotonic_checks':0,'cycle_edges':0}
    def primal(cert,ids,t,context_specific=False):
        rr=[rows[i] for i in ids]
        xx=[x[i%66] if not context_specific else (x[i%66]+[0]*92 if i<66 else [0]*92+x[i%66]) for i in ids]
        intervals=[(lu[z['id'],z['cell']][0],t*lu[z['id'],z['cell']][1]) for z in rr]
        forced=[[i,j] for i in range(len(rr)) for j in range(len(rr)) if rr[i]['cell']==rr[j]['cell'] and intervals[j][1]<intervals[i][0]]
        assert forced==cert['forced_pairs'] and cert['status']=='feasible'
        beta=list(map(Q,cert['beta'])); assert len(beta)==len(xx[0])
        zz=[sum(a*b for a,b in zip(row,beta)) for row in xx]
        assert zz==list(map(Q,cert['latent']))
        margins=[zz[i]-zz[j] for i,j in forced]
        assert margins==list(map(Q,cert['margins'])) and all(v>=1 for v in margins)
        yy=list(map(Q,cert['responses']));assert len(yy)==len(rr)
        assert all(0<v and lo<=v<=hi for v,(lo,hi) in zip(yy,intervals))
        count=0
        for i in range(len(rr)):
            for j in range(len(rr)):
                if rr[i]['cell']==rr[j]['cell'] and zz[i]<=zz[j]:
                    assert yy[i]<=yy[j];count+=1
        totals['primal_cases']+=1;totals['forced_comparisons']+=len(forced)
        totals['responses']+=len(yy);totals['isotonic_checks']+=count
        return xx
    conclusions=[]
    for e in r['exons']:
        ids=[i for i,z in enumerate(rows) if z['exon']==e['exon']]
        assert ids==e['row_indices'];t=Q(e['M0_t']);assert t==Q(e['arbitrary_order_t']) and t>=1
        xx=primal(e['at_arbitrary_boundary'],ids,t)
        if t>1:
            cycle=e['lower_cycle'];assert cycle
            weights=[]
            for k,edge in enumerate(cycle):
                assert edge['after']==cycle[(k+1)%len(cycle)]['before']
                ratios={c:lu[edge['after'],c][0]/lu[edge['before'],c][1] for c in ('K562','KBM7')}
                w=max(ratios.values());assert w==Q(edge['weight'])
                assert sorted(c for c,v in ratios.items() if v==w)==sorted(edge['contexts'])
                weights.append(w);totals['cycle_edges']+=1
            assert min(weights)==t
        else: assert not e['lower_cycle']
        one=[row for row,i in zip(xx,ids) if i<66];rank=Matrix(one).rank()
        if e['exon']==27: assert rank==10==e['rank']
        conclusions.append({'exon':e['exon'],'t':str(t),'rank':rank,'local_states':len(one),'same_as_arbitrary_order':True})
    for e in r['m1_at_one']:
        ids=[i for i,z in enumerate(rows) if z['exon']==e['exon']]
        assert ids==e['row_indices'];primal(e['certificate'],ids,Q(1),True)
    assert len(conclusions)==5 and totals['primal_cases']==10
    assert all(sha(files[k])==v for k,v in pins.items())
    out=Path(__file__).with_name('verification.json')
    if out.exists(): raise FileExistsError('Preserve completed review')
    report={'status':'PASS_INDEPENDENT_EXACT_BOUNDARY_CERTIFICATES','pins':pins,'counts':totals,'conclusions':conclusions,'checker_sha256':sha(Path(__file__)),'scope':'Separate matrix reconstruction and rational primal/cycle check without optimization. Numerical verification, not independent biological replication or calibrated uncertainty.'}
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(totals))
if __name__=='__main__':main()
