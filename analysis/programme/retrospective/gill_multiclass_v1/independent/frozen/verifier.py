"""Independent Gill verifier. Preparation reads metadata only.

The evaluate phase must not be invoked until root authorizes validation access.
No author scientific implementation is imported in either phase.
"""
from pathlib import Path
import argparse
import collections
from datetime import datetime, timezone
import gzip
import hashlib
import itertools
import json
import math
import platform
import sys
import time
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT/'outputs/synonymous-rna-accessibility/analysis/programme/retrospective/gill_multiclass_v1'
TRAIN = ROOT/'work/gill_multiclass_cycle3/training_v1'
EXPECTED_PROTOCOL = '4ab4e85c14203c0ab8c983e8433c766f1566f373a335fd3fda4fb3981371c1aa'
EXPECTED_INVENTORY = '2e00b698dee109e357dc01435d5eacf8b341388f36463970454ade230664f350'
EXPECTED_TRAIN_RECEIPT = 'c6fb3309a93635efab38f1469516617c3f1b4f857c88b4fbfff91be21115bdef'
SCENARIOS = list(itertools.product(range(3),repeat=2))
PAIR45 = list(itertools.combinations_with_replacement(range(9),2))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_json(p):
    raw=Path(p).read_bytes()
    return json.loads(gzip.decompress(raw) if str(p).endswith('.gz') else raw)
def dump(p,x):
    raw=(json.dumps(x,sort_keys=True,allow_nan=False,separators=(',',':'))+'\n').encode()
    Path(p).write_bytes(gzip.compress(raw,mtime=0) if str(p).endswith('.gz') else raw)
def cosine(a,b):
    aa=float(np.dot(a,a));bb=float(np.dot(b,b))
    return None if aa==0 or bb==0 else float(np.dot(a,b)/math.sqrt(aa*bb))
def codons(key):return (key[:1],key[1:4],key[4:7],key[7:8])

def setup_class(c):
    keys=c['observed_keys'];n=len(keys);labels=[codons(k) for k in keys]
    assert len(set(keys))==n and all(len(k)==8 and set(k)<=set('ACGT') for k in keys)
    levels=[sorted({t[j] for t in labels}) for j in range(4)]
    assert set(labels)==set(itertools.product(*levels))
    assert n==c['expected_vertices']==c['observed_vertices'] and c['complete']
    X=np.ones((n,1+sum(len(l)-1 for l in levels)))
    col=1
    for j,L in enumerate(levels):
        for level in L[1:]:X[:,col]=[t[j]==level for t in labels];col+=1
    rank=int(np.linalg.matrix_rank(X))
    assert rank==X.shape[1]==c['codon_main_effect_rank']
    assert n-rank==c['codon_residual_dimensions']>=3
    assert list(map(len,levels))==c['codon_factor_sizes']
    adjacency=np.zeros((n,n),dtype=bool)
    for i,k in enumerate(keys):
        for j in range(i):
            if sum(a!=b for a,b in zip(k,keys[j]))==1:adjacency[i,j]=adjacency[j,i]=True
    assert int(adjacency.sum())//2==c['expected_edges']==c['observed_edges']
    reached={0};todo=[0]
    while todo:
        v=todo.pop()
        for w in np.flatnonzero(adjacency[v]):
            if int(w) not in reached:reached.add(int(w));todo.append(int(w))
    assert len(reached)==n and c['expected_component_sizes']==[n]
    assert not c['reference_protein'] and c['central_amino_acid_stratum_157_158']!='VP'
    return keys,X,adjacency

def project(X,Y):
    # Fit an indicator matrix, rather than implementing conditional means.
    if Y.ndim==1:Y=Y[:,None]
    centered=Y-Y.mean(axis=0)
    fitted=X@np.linalg.lstsq(X,centered,rcond=None)[0]
    Q=centered-fitted
    return Y-Q,Q

def floyd_barriers(adjacency,F):
    """All pair capacities for every row field F by max-min matrix closure."""
    m,n=F.shape
    endmins=np.minimum(F[:,:,None],F[:,None,:])
    allowed=adjacency|np.eye(n,dtype=bool)
    C=np.where(allowed[None,:,:],endmins,-np.inf)
    for k in range(n):
        C=np.maximum(C,np.minimum(C[:,:,k,None],C[:,None,k,:]))
    assert np.isfinite(C).all()
    B=endmins-C
    return B.max(axis=(1,2)),[[int(v) for v in np.unravel_index(int(np.argmax(b)),b.shape)] for b in B]

def midranks(x):
    order=np.argsort(x,kind='stable');r=np.empty(len(x),dtype=float);i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]:j+=1
        r[order[i:j]]=(i+j-1)/2+1;i=j
    return r-r.mean()

def get_metadata():
    assert sha(BASE/'protocol.json')==EXPECTED_PROTOCOL
    raw=gzip.decompress((BASE/'inventory.json.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest()==EXPECTED_INVENTORY
    p=read_json(BASE/'protocol.json');inv=json.loads(raw)
    ids={part:set(p['class_ids_by_partition'][part]) for part in ['development','validation']}
    assert len(ids['development'])==1124 and len(ids['validation'])==1052
    assert not (ids['development']&ids['validation'])
    classes=[c for c in inv['classes'] if c['class_id'] in ids['development']|ids['validation']]
    assert len(classes)==2176 and len({c['class_id'] for c in classes})==2176
    flatkeys=[k for c in classes for k in c['observed_keys']]
    assert len(flatkeys)==len(set(flatkeys))==38338
    strata={part:{c['central_amino_acid_stratum_157_158'] for c in classes if c['class_id'] in ids[part]} for part in ids}
    assert not strata['development']&strata['validation']
    # Assigned possible strata are200/199; eligibility leaves162/154 represented.
    assert len(strata['development'])==162 and len(strata['validation'])==154
    for part in ids:
        assert all(c['candidate_partition']==part for c in classes if c['class_id'] in ids[part])
        expected=p['development_vertices'] if part=='development' else p['primary_vertices']
        assert sum(c['observed_vertices'] for c in classes if c['class_id'] in ids[part])==expected
    assert p['pseudocounts']==[0.25,0.5,1.0]
    return p,classes,ids

def self_test():
    p,classes,ids=get_metadata();architectures=collections.Counter();max_error=0.;synthetic_rows=[]
    for c in classes:
        keys,X,A=setup_class(c);architectures[tuple(c['codon_factor_sizes'])]+=1
        labels=[codons(k) for k in keys]
        # Abstract deterministic numeric fields only, no observed counts.
        Y=np.column_stack([np.array([((i+1)*(j+2))%11 for i in range(len(keys))],dtype=float) for j in range(6)])
        P,Q=project(X,Y)
        means=Y.mean(axis=0)
        conditional=sum(np.vstack([Y[[t[j]==u[j] for t in labels]].mean(axis=0) for u in labels]) for j in range(4))-3*means
        error=float(np.max(np.abs(P-conditional)));max_error=max(max_error,error)
        assert error<1e-10
        assert np.max(np.abs(X.T@Q))<1e-9
        if architectures[tuple(c['codon_factor_sizes'])]==1:
            synthetic_rows.append(class_metrics(c,X,A,Y,np.exp2(Y[:,:3]),0.5))
    # Hand-check a max-min path and a width-one projection counterexample.
    A=np.array([[False,True,False],[True,False,True],[False,True,False]])
    b,_=floyd_barriers(A,np.array([[2.,0.,2.],[1.,1.,1.]]));assert np.array_equal(b,[2.,0.])
    Z=np.arange(36,dtype=float).reshape(6,6)/7;G=Z@Z.T
    for u,v in PAIR45:
        i,j=SCENARIOS[u];k,l=SCENARIOS[v]
        terms=[G[i,k],G[3+j,3+l],-G[i,3+l],-G[3+j,k]]
        assert abs(sum(terms)-np.dot(Z[i]-Z[3+j],Z[k]-Z[3+l]))<1e-10
    assert len(PAIR45)==45
    synthetic_aggregate=pooled(synthetic_rows)
    json.dumps(synthetic_aggregate,allow_nan=False)
    return {'status':'metadata_and_synthetic_checks_passed','classes':len(classes),'vertices':sum(c['observed_vertices'] for c in classes),
      'architectures':{str(k):v for k,v in architectures.items()},'largest_class':max(c['observed_vertices'] for c in classes),
      'sum_n_cubed':sum(c['observed_vertices']**3 for c in classes),'maximum_indicator_vs_conditional_projection_error':max_error,
      'synthetic_all_metric_architectures':len(synthetic_rows),'actual_counts_opened':False,'validation_artifacts_opened':False}

def class_metrics(c,X,adj,L,rawRNA,cp):
    P,Q=project(X,L)
    R,D=Q[:,:3],Q[:,3:]
    gram=Q.T@Q
    T=(L[:,0]+L[:,1]-L[:,3]-L[:,4])/2;V=L[:,2]-L[:,5]
    qt=(Q[:,0]+Q[:,1]-Q[:,3]-Q[:,4])/2;qv=Q[:,2]-Q[:,5]
    pt=T-qt;pv=V-qv;tc=T-T.mean();vc=V-V.mean()
    scenarios=np.column_stack([L[:,i]-L[:,3+j] for i,j in SCENARIOS])
    qs=np.column_stack([Q[:,i]-Q[:,3+j] for i,j in SCENARIOS]);ps=scenarios-qs
    terms45=[]
    for u,v in PAIR45:
        i,j=SCENARIOS[u];k,l=SCENARIOS[v]
        terms=[float(gram[i,k]),float(gram[3+j,3+l]),float(-gram[i,3+l]),float(-gram[3+j,k])]
        direct=float(np.dot(qs[:,u],qs[:,v]));delta=direct-math.fsum(terms)
        bound=2e-10*max(sum(abs(x) for x in terms),np.linalg.norm(Q)**2,np.finfo(float).tiny)
        assert abs(delta)<=bound
        terms45.append({'scenarios':[u,v],'terms':terms,'total':direct,'identity_residual':delta})
    rb=(R[:,0]+R[:,1])/2;db=(D[:,0]+D[:,1])/2
    primary_terms=[float(np.dot(rb,R[:,2])),float(np.dot(db,D[:,2])),float(-np.dot(rb,D[:,2])),float(-np.dot(db,R[:,2]))]
    assert abs(math.fsum(primary_terms)-np.dot(qt,qv))<=2e-10*max(np.linalg.norm(Q)**2,np.finfo(float).tiny)
    fields=np.vstack([scenarios.T,ps.T,T[None,:],pt[None,:]])
    barriers,argmax=floyd_barriers(adj,fields)
    stats={'residual_TT':float(np.dot(qt,qt)),'residual_TV':float(np.dot(qt,qv)),'residual_VV':float(np.dot(qv,qv)),
      'residual_cosine':cosine(qt,qv),'residual_mse':float(np.mean((qt-qv)**2)),
      'zero_residual_mse':float(np.mean(qv**2)),'whole_TT':float(np.dot(tc,tc)),
      'whole_TV':float(np.dot(tc,vc)),'whole_VV':float(np.dot(vc,vc)),
      'whole_cosine':cosine(tc,vc),'whole_mse':float(np.mean((tc-vc)**2)),
      'codon_main_prediction_mse':float(np.mean((pt-pt.mean()-vc)**2)),
      'constant_prediction_mse':float(np.mean(vc**2)),'projection_rank':int(X.shape[1]),'residual_dimension':len(T)-X.shape[1]}
    pair_stats={}
    for i,j in itertools.combinations(range(3),2):
        pair_stats[str((i,j))]={'RNA_residual_cosine':cosine(R[:,i],R[:,j]),'DNA_residual_cosine':cosine(D[:,i],D[:,j]),
          'RNA_midrank_correlation':cosine(midranks(rawRNA[:,i]),midranks(rawRNA[:,j]))}
    return {'class_id':c['class_id'],'n':len(T),'stratum':c['central_amino_acid_stratum_157_158'],
      'all_factor_cliques':c['codon_factors_all_cliques'],'stats':stats,'primary_terms':primary_terms,
      'decomposition45':terms45,'RNA_residual_gram':gram[:3,:3].tolist(),'DNA_residual_gram':gram[3:,3:].tolist(),
      'RNA_DNA_residual_gram':gram[:3,3:].tolist(),'pair_stats':pair_stats,'scenario_barriers':barriers[:9].tolist(),
      'scenario_codon_barriers':barriers[9:18].tolist(),'training_barrier':float(barriers[18]),
      'training_codon_barrier':float(barriers[19]),'heldout_barrier':float(barriers[8]),
      'heldout_codon_barrier':float(barriers[17]),'barrier_argmax_vertex_indices':argmax,
      'independent_training_fields':{'T':T.tolist(),'PT':pt.tolist(),'QT':qt.tolist()},
      'raw_log_input_norm':float(np.linalg.norm(L)),
      'centered_log_input_norm':float(np.linalg.norm(L-L.mean(axis=0)))}

def pooled(rows):
    terms={k:math.fsum(r['stats'][k]/r['n'] for r in rows) for k in ['residual_TT','residual_TV','residual_VV','whole_TT','whole_TV','whole_VV']}
    def rho(t,prefix):
        a,b=t[prefix+'TT'],t[prefix+'VV']
        return None if a==0 or b==0 else t[prefix+'TV']/math.sqrt(a*b)
    loss={k:math.fsum(r['stats'][k] for r in rows) for k in ['residual_mse','zero_residual_mse','whole_mse','codon_main_prediction_mse','constant_prediction_mse']}
    group=collections.defaultdict(list)
    for r in rows:group[r['stratum']].append(r)
    influence=[]
    for s,rr in sorted(group.items()):
        left={k:terms[k]-math.fsum(r['stats'][k]/r['n'] for r in rr) for k in terms}
        influence.append({'stratum':s,'class_count':len(rr),'pooled_residual_cosine_without':rho(left,'residual_')})
    energy_share={k:None if terms[k]==0 else max(math.fsum(r['stats'][k]/r['n'] for r in rr)/terms[k] for rr in group.values()) for k in ['residual_TT','residual_VV']}
    pairs45=[]
    for i in range(45):
        pairs45.append({'scenarios':PAIR45[i],'terms':[math.fsum(r['decomposition45'][i]['terms'][j]/r['n'] for r in rows) for j in range(4)],
          'total':math.fsum(r['decomposition45'][i]['total']/r['n'] for r in rows)})
    primary=[math.fsum(r['primary_terms'][j]/r['n'] for r in rows) for j in range(4)]
    diag={p['scenarios'][0]:p['total'] for p in pairs45 if p['scenarios'][0]==p['scenarios'][1]}
    for p in pairs45:
        a,b=p['scenarios'];p['cosine']=None if diag[a]==0 or diag[b]==0 else p['total']/math.sqrt(diag[a]*diag[b])
    grams={kind:[[math.fsum(r[kind][i][j]/r['n'] for r in rows) for j in range(3)] for i in range(3)]
      for kind in ['RNA_residual_gram','DNA_residual_gram','RNA_DNA_residual_gram']}
    pair_cosines={kind:{str((i,j)):None if g[i][i]==0 or g[j][j]==0 else g[i][j]/math.sqrt(g[i][i]*g[j][j])
      for i,j in itertools.combinations(range(3),2)} for kind,g in grams.items() if kind!='RNA_DNA_residual_gram'}
    influence_values=[x['pooled_residual_cosine_without'] for x in influence if x['pooled_residual_cosine_without'] is not None]
    return {'classes':len(rows),'vertices':sum(r['n'] for r in rows),'strata':len(group),
      'terms':terms,'pooled_residual_cosine':rho(terms,'residual_'),'pooled_whole_cosine':rho(terms,'whole_'),
      'operational_R_ge_0_5':None if rho(terms,'residual_') is None else rho(terms,'residual_')>=0.5,
      'sum_class_mean_losses':loss,'primary_decomposition_terms':primary,'decomposition45':pairs45,
      'leave_one_stratum_out':influence,'leave_one_stratum_out_range':None if not influence_values else [min(influence_values),max(influence_values)],
      'maximum_stratum_energy_fraction':energy_share,'residual_grams':grams,'RNA_DNA_pooled_pair_cosines':pair_cosines,
      'undefined_class_residual_cosines':sum(r['stats']['residual_cosine'] is None for r in rows),
      'residual_loss_gain_vs_zero':None if loss['zero_residual_mse']==0 else 1-loss['residual_mse']/loss['zero_residual_mse'],
      'barrier_extrema_by_clique_status':{str(flag):{kind:max((r[kind] for r in rows if r['all_factor_cliques']==flag),default=None)
        for kind in ['training_barrier','training_codon_barrier','heldout_barrier','heldout_codon_barrier']} for flag in [False,True]}}

def evaluate(args):
    # This phase is explicitly gated by root's future validation receipt/hash.
    frozen=read_json(HERE/'freeze.json');assert sha(__file__)==frozen['verifier_sha256']
    assert args.validation_receipt and args.validation_receipt_sha and args.validation_counts and args.validation_counts_sha
    assert sha(args.validation_receipt)==args.validation_receipt_sha
    assert sha(args.validation_counts)==args.validation_counts_sha
    assert sha(TRAIN/'receipt.json')==EXPECTED_TRAIN_RECEIPT
    receipt=read_json(TRAIN/'receipt.json');assert not receipt['D_opened']
    assert receipt['inputs']['protocol_sha256']==EXPECTED_PROTOCOL
    assert receipt['inputs']['inventory_sha256']==EXPECTED_INVENTORY
    for name in ['counts.json.gz','predictions.json.gz']:
        assert sha(TRAIN/name)==receipt['outputs'][name]
    p,classes,ids=get_metadata();train=read_json(TRAIN/'counts.json.gz');held=read_json(args.validation_counts)
    keys={k for c in classes for k in c['observed_keys']}
    assert set(train)==keys==set(held)
    for k in keys:
        assert set(train[k])==set('BCEFG') and set(held[k])==set('DG')
        assert train[k]['G']==held[k]['G']
        assert all(type(v) is int and v>=0 for v in [*train[k].values(),*held[k].values()])
        assert sum(train[k][j] for j in 'EFG')>50
    result={};start=time.perf_counter()
    for cp in p['pseudocounts']:
        parts={part:[] for part in ids}
        for c in classes:
            kk,X,A=setup_class(c)
            counts=np.array([[train[k]['B'],train[k]['C'],held[k]['D'],train[k]['E'],train[k]['F'],train[k]['G']] for k in kk],dtype=float)
            L=np.log2(counts+cp);assert np.isfinite(L).all()
            row=class_metrics(c,X,A,L,counts[:,:3],cp)
            part='validation' if c['class_id'] in ids['validation'] else 'development';parts[part].append(row)
        result[str(cp)]={'classes':parts,'aggregate':{part:pooled(rr) for part,rr in parts.items()}}
    assert set(result)==set(['0.25','0.5','1.0'])
    output=Path(args.output);assert not output.exists();output.mkdir(parents=True)
    dump(output/'independent_results.json.gz',result)
    dump(output/'receipt.json',{'status':'computed_independently_no_author_metrics_read','runtime_seconds':time.perf_counter()-start,
      'protocol_sha256':EXPECTED_PROTOCOL,'inventory_sha256':EXPECTED_INVENTORY,'training_receipt_sha256':EXPECTED_TRAIN_RECEIPT,
      'validation_receipt_sha256':sha(args.validation_receipt),'validation_counts_sha256':sha(args.validation_counts),
      'source_sha256':sha(__file__),'result_sha256':sha(output/'independent_results.json.gz'),
      'python':platform.python_version(),'numpy':np.__version__,'all_classes':len(classes),'all_pseudocounts':p['pseudocounts']})
    print('Independent metrics saved; author comparison requires a separate explicit read.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['prepare','evaluate'],required=True)
    ap.add_argument('--validation-receipt');ap.add_argument('--validation-receipt-sha');ap.add_argument('--validation-counts');ap.add_argument('--validation-counts-sha');ap.add_argument('--output')
    args=ap.parse_args()
    if args.phase=='prepare':
        target=Path(args.output) if args.output else HERE/'preparation_checks.json'
        assert not target.exists();dump(target,self_test());print('Metadata/synthetic preparation passed; no observed counts opened.')
    else:evaluate(args)

if __name__=='__main__':main()
