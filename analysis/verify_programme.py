"""Independent, no-folding replay of the discovery audit's consequential outputs."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
import hashlib
import json
import math
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/exact'))
from independent_oracle import CODONS
TRANSLATE={c:a for a,cc in CODONS.items() for c in cc}
POINTS=list(product((0,1),repeat=3))


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def graph(strings):
    return [[j for j,b in enumerate(strings) if sum(x!=y for x,y in zip(a,b))==1] for a in strings]


def components(g,active):
    labels=[-1]*len(g);count=0
    for i in active:
        if labels[i]>=0:continue
        stack=[i];labels[i]=count
        while stack:
            for j in g[stack.pop()]:
                if j in active and labels[j]<0:labels[j]=count;stack.append(j)
        count+=1
    return labels


def endpoint_B(g,q,a,b):
    for cutoff in sorted(set(q),reverse=True):
        if min(q[a],q[b])<cutoff:continue
        cc=components(g,{i for i,v in enumerate(q) if v>=cutoff})
        if cc[a]==cc[b]:return min(q[a],q[b])-cutoff
    raise AssertionError('Disconnected baseline endpoints')


def pair_list(dbn):
    stack=[];pairs=[]
    for i,c in enumerate(dbn):
        if c=='(':stack.append(i)
        elif c==')':pairs.append((stack.pop(),i))
    assert not stack
    return pairs


def check_manifests():
    count=0
    for directory in ['design_gate_v1','overlap_probe_v1','overlap_positive_v1','robust_certificates_v1']:
        data=json.loads((ROOT/'results/programme'/directory/'manifest.json').read_text())
        assert data['run']['status']=='completed'
        for row in data['input_artifacts']:
            assert sha(ROOT/row['path'])==row['sha256']==row['sha256_after']
            count+=1
        for row in data['outputs']:
            assert sha(ROOT/row['path'])==row['sha256'];count+=1
    return count


def check_predictions():
    path=ROOT/'results/programme/cube_pilot_v1/predictions.json'
    data=json.loads(path.read_text());feas=json.loads(path.with_name('feasibility.json').read_text())
    assert feas['prediction_sha256']==sha(path)
    assert data['partition_functions_run']==0
    assert not path.with_name('development.json').exists()
    assert not path.with_name('heldout.json').exists()
    for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest
    g=graph(POINTS);checked=nonzero=unate=0;max_B=0.;max_theta=0.
    for scaffold in data['scaffolds']:
        for variant,pred in scaffold['predictions'].items():
            sequences=pred['sequences']
            assert graph(sequences)==g
            assert len({''.join(TRANSLATE[s[i:i+3]] for i in range(0,len(s),3)) for s in sequences})==1
            assert all(s[i]+s[j] in {'AU','UA','GC','CG','GU','UG'} for s in sequences for i,j in scaffold['pairs0'])
            pairsets=[pair_list(x) for x in pred['structure_library']]
            target=set(map(tuple,scaffold['pairs0']))
            occupancy=[len(set(pp)&target)/9 for pp in pairsets]
            q=[]
            for seq,energy in zip(sequences,pred['library_energies_kcal']):
                for pp,e in zip(pairsets,energy):
                    assert (e is not None)==all(seq[i]+seq[j] in {'AU','UA','GC','CG','GU','UG'} for i,j in pp)
                # Separate log-sum-exp implementation, from frozen energy tables.
                weights=[-e/(.00198720425864083*310.15) if e is not None else -math.inf for e in energy]
                m=max(weights);normalizer=sum(math.exp(w-m) for w in weights)
                q.append(sum(math.exp(w-m)*p for w,p in zip(weights,occupancy))/normalizer)
            f=pred['finite_competitor']
            assert max(abs(a-b) for a,b in zip(q,f['scores']))<1e-12
            assert abs(endpoint_B(g,q,0,6)-f['B_primary'])<1e-12
            for row in f['all_endpoint_lifetimes']:
                assert abs(endpoint_B(g,q,row['a'],row['b'])-row['B'])<1e-12
            theta=sum((-1)**(3-sum(x))*v for x,v in zip(POINTS,q))
            assert abs(theta-f['mobius']['111'])<1e-12
            nonzero+=abs(theta)>1e-10;max_theta=max(max_theta,abs(theta))
            max_B=max(max_B,max(endpoint_B(g,q,a,b) for a,b in combinations(range(8),2)))
            is_unate=all(min(d:=[q[i+(1<<(2-k))]-q[i] for i,x in enumerate(POINTS) if not x[k]])>=-1e-12 or max(d)<=1e-12 for k in range(3))
            unate+=is_unate;checked+=1
        assert scaffold['maximum100nt_identity']<.8 and scaffold['maximum_old90nt_identity']<.8
    # Global epistasis control: nonlinear transformation alone creates theta123.
    nonlinear=[F(sum(x)**3,27) for x in POINTS]
    assert sum((-1)**(3-sum(x))*v for x,v in zip(POINTS,nonlinear))==F(2,9)
    assert all(endpoint_B(g,nonlinear,a,b)==0 for a,b in combinations(range(8),2))
    assert feas['summary']['development']['covered_intervention_contrasts']==0
    return {'landscapes':checked,'score_reconstruction_max_tolerance':1e-12,
            'nonzero_three_way_coefficients_gt1e_10':nonzero,'coordinate_unate_fields_tolerance1e_12':unate,
            'maximum_absolute_three_way':max_theta,'maximum_any_endpoint_lifetime':max_B,
            'partition_functions_replayed':0,'heldout_outcomes_unopened':True,
            'global_epistasis_counterexample_theta123':'2/9'}


def check_robust_boxes():
    data=json.loads((ROOT/'results/programme/robust_certificates_v1/results.json').read_text())
    positives=[]
    for row in data['rows']:
        name=row['id']+'_'+row['condition']+'.json'
        rs=json.loads((ROOT/'results/expansion/rnastructure'/name).read_text())
        vn=json.loads((ROOT/'results/native'/name).read_text())
        assert rs['genotypes']==vn['genotypes']
        g=graph(vn['genotypes']);base=components(g,set(range(len(g))))
        lower=[min(F(str(v)),F(str(bounds[0]))) for v,bounds in zip(vn['pair_retention'],rs['retention_rounding_bounds'])]
        upper=[max(F(str(v)),F(str(bounds[1]))) for v,bounds in zip(vn['pair_retention'],rs['retention_rounding_bounds'])]
        positive=False
        values=sorted(set(lower+upper+[F(0),F(1)]))
        for lo,hi in zip(values,values[1:]):
            t=(lo+hi)/2;possible=components(g,{i for i,v in enumerate(upper) if v>=t})
            seen={}
            for i,v in enumerate(lower):
                if v<t:continue
                if base[i] in seen and seen[base[i]]!=possible[i]:positive=True;break
                seen[base[i]]=possible[i]
            if positive:break
        assert positive==row['any_certified_positive_barrier']
        if positive:
            w=row['witness'];a,b=w['a'],w['b']
            lo,hi=map(F,w['rational_barrier_interval_open_closed']);t=(lo+hi)/2
            possible=components(g,{i for i,v in enumerate(upper) if v>=t})
            assert base[a]==base[b] and lower[a]>=hi and lower[b]>=hi
            assert hi>lo and possible[a]!=possible[b]
            positives.append((row['id'],row['condition']))
    return {'graphs':len(data['rows']),'positive_condition_rows':len(positives),
            'distinct_positive_blocks':len({x[0] for x in positives})}


if __name__=='__main__':
    result={'manifest_artifact_hash_checks':check_manifests(),
            'cube_predictions':check_predictions(),'robust_certificates':check_robust_boxes()}
    print(json.dumps(result,indent=2))
