"""Separately authored certificate replay, without fitting or importing D code."""
import argparse, json, hashlib
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from sympy import Matrix

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if not __debug__:raise RuntimeError('Assertions required')
    p=argparse.ArgumentParser();p.add_argument('--design',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('No overwrite')
    files={'map':a.design/'exact_v1/mapping.json','matrix':a.design/'inventory_v1/nested_matrices.json','certificates':a.design/'exact_v1/certificates.json','result':a.design/'nested_model/actual_v1/result.json','source':a.source}
    pins={'map':'a2f428f5308cb78a827bfbd23f7d7756b98980e3585f0ccf6a733859c2388436','matrix':'665c4d0586310d09868b07f01163875e82708c3441201ab7757d7c897b3ab874','certificates':'64297d0dcc537a1fd7304af4a925fcca8c01e6a9ac2e9c39af1919ba9a526534','result':'9b9726f456a0309412df34530ec03b826e70f7ac5610617142ac96cf6e827e2e','source':'eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'}
    assert all(sha(files[k])==v for k,v in pins.items())
    objects={k:json.loads(f.read_text()) for k,f in files.items()};mapping=objects['map'];models=objects['matrix'];result=objects['result'];source=objects['source'];old=objects['certificates']
    assert len(mapping)==66 and len({r['id'] for r in mapping})==66
    mapping.sort(key=lambda r:r['id']);exons=sorted({r['exon'] for r in mapping});alleles=sorted({al for r in mapping for al in r['alleles']})
    cols=[f'exon:{e}' for e in exons]+alleles
    assert cols==old['column_labels']
    x=[[int(r['exon']==e) for e in exons]+[int(al in r['alleles']) for al in alleles] for r in mapping]
    assert x==old['matrix'] and [r['id'] for r in mapping]==old['row_ids']
    # Authenticate that the supplied counted words encode exactly the two
    # descriptors. No reference imputation outside these words is made.
    for r in mapping:
        observed=r['observed_DNA_word'];reference=r['reference_word'];start=r['observed_start']
        assert len(observed)==len(reference)==r['observed_end']-start+1
        diffs={(start+k,b,c) for k,(b,c) in enumerate(zip(reference,observed)) if b!=c}
        assert diffs=={(d['pos'],d['ref'],d['alt']) for d in r['changes']}
        assert r['observed_bases']=={str(start+k):b for k,b in enumerate(observed)}
    X=Matrix(x);rank=X.rank();left=X.T.nullspace();assert rank==65 and len(left)==1
    support=[mapping[i]['id'] for i,z in enumerate(left[0]) if z]
    ranks={str(e):Matrix([row for row,r in zip(x,mapping) if r['exon']==e]).rank() for e in exons}
    allrows=[{'cell':c,'id':r['id'],'exon':r['exon']} for c in ['K562','KBM7'] for r in mapping]
    assert allrows==models['rows']
    m0=x+x;m1=[r+[0]*92 for r in x]+[[0]*92+r for r in x]
    xi=[r+[int(i==j) for j in range(66)] for i,r in enumerate(x)]
    m2=[r+[0]*158 for r in xi]+[[0]*158+r for r in xi]
    expected={'strict_shared':m0,'context_specific_additive':m1,'context_specific_all_observed_pairs':m2}
    assert expected==models['matrices']
    intervals={(r['state'],r['cell']):(Q(r['lower']),None if r['upper'] is None else Q(r['upper'])) for r in source['state_cells']}
    assert len(intervals)==132
    expected_cases={(str(g),e,m,arm) for g in [Q(1),Q(5,4),Q(3,2),Q(2)] for e in exons for m in expected for arm in (['both'] if m=='strict_shared' else ['both','K562','KBM7'])}
    seen=set();checked_pairs=checked_responses=checked_dual_rows=0
    classifications=[]
    for r in result['results']:
        key=(r['gamma'],r['exon'],r['model'],r['arm']);assert key not in seen;seen.add(key)
        ids=[i for i,z in enumerate(allrows) if z['exon']==r['exon'] and (r['arm']=='both' or z['cell']==r['arm'])]
        assert ids==r['row_indices'];rows=[allrows[i] for i in ids];design=[expected[r['model']][i] for i in ids];g=Q(r['gamma'])
        lu=[intervals[(z['id'],z['cell'])] for z in rows];lower=[l/g for l,u in lu];upper=[None if u is None else u*g for l,u in lu]
        pairs=[[i,j] for i in range(len(rows)) for j in range(len(rows)) if rows[i]['cell']==rows[j]['cell'] and upper[i] is not None and upper[i]<lower[j]]
        assert pairs==r['forced_pairs'] and len(pairs)==r['forced_pair_count'];checked_pairs+=len(pairs)
        delta=[[b-a for a,b in zip(design[i],design[j])] for i,j in pairs]
        if r['status']=='feasible':
            beta=list(map(Q,r['beta']));latent=[sum(Q(v)*w for v,w in zip(row,beta)) for row in design]
            assert latent==list(map(Q,r['latent']))
            gaps=[sum(Q(v)*w for v,w in zip(row,beta)) for row in delta]
            assert gaps==list(map(Q,r['margins'])) and all(z>=1 for z in gaps)
            y=list(map(Q,r['responses']))
            assert len(y)==len(rows) and all(y[i]>0 and lower[i]<=y[i] and (upper[i] is None or y[i]<=upper[i]) for i in range(len(y)))
            assert all(y[i]<=y[j] for i in range(len(y)) for j in range(len(y)) if rows[i]['cell']==rows[j]['cell'] and latent[i]<=latent[j])
            checked_responses+=len(y)
        elif r['status']=='infeasible':
            ww=[(z['constraint_index'],Q(z['weight'])) for z in r['dual']]
            assert sum(w for k,w in ww)==1 and all(w>0 for k,w in ww)
            assert len({k for k,w in ww})==len(ww)
            assert all(z['difference']==delta[z['constraint_index']] for z in r['dual'])
            assert all(sum(w*delta[k][j] for k,w in ww)==0 for j in range(len(design[0])))
            checked_dual_rows+=len(ww)
        else:raise ValueError('Unresolved result')
        classifications.append({'gamma':r['gamma'],'exon':r['exon'],'model':r['model'],'arm':r['arm'],'status':r['status']})
    assert seen==expected_cases and len(seen)==140
    assert all(sha(files[k])==v for k,v in pins.items())
    report={'status':'PASS_SEPARATE_EXACT_CERTIFICATE_AND_DESIGN_RECONSTRUCTION','input_hashes':pins,'rank':rank,'per_exon_rank':ranks,'left_null_support_ids':support,'parameter_nullity':len(cols)-rank,'matrix_entries_checked':sum(len(r) for x0 in expected.values() for r in x0),'cases':len(seen),'forced_pairs_checked':checked_pairs,'feasible_response_values_checked':checked_responses,'dual_rows_checked':checked_dual_rows,'status_counts':dict(Counter(r['status'] for r in classifications)),'infeasible_cases':[r for r in classifications if r['status']=='infeasible'],'all_classifications':classifications,'checker_sha256':sha(Path(__file__)),'limits':'Checks pinned local-word mapping and certificates with independent matrix construction and SymPy rank, without refitting. Not independent source extraction, biological replication, calibration or full-haplotype authentication.'}
    a.out.mkdir(parents=True);(a.out/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:report[k] for k in ['status','rank','per_exon_rank','cases','forced_pairs_checked','feasible_response_values_checked','dual_rows_checked','status_counts']}))

if __name__=='__main__':main()
