"""Independent exact observed robustness verification; no production imports or LP."""
import sys
if not __debug__:raise RuntimeError('Optimized Python is not allowed')
from pathlib import Path
from fractions import Fraction as Q
import argparse
import datetime
import gzip
import hashlib
import json
import platform
import time

RESULT_SHA='861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa'
PRIMARY_SHA='3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466'
GZIP_SHA='e3c3d061052643046de58dc789866899b0fb14ade83f06170e2cf32e9e54e145'
FREEZE_SHA='340624e795e600abbbaeb3daefa89b054ceca7e356855014d320df20b043a10f'
CONTRACT_SHA='f0fb3d2ced5dc36d0794eb76b7a257351534f8dd3a784f33c024be4b500a9c38'
WT='GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'

def need(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def read_pinned(path,digest):
    raw=path.read_bytes();need(sha(raw)==digest,'Hash mismatch: '+str(path));return raw

def design(sequences):
    states=[[s[k:k+3] for k in (23,38,29,32,41)] for s in sequences]
    levels=[sorted({row[j] for row in states}) for j in range(5)]
    X=[[1,*[int(row[j]==level) for j in range(5) for level in levels[j][1:]]] for row in states]
    need(all(len(row)==8 for row in X),'Feature design size')
    return X,states

def certify(cert,X,L,U):
    n=len(X);need(n==len(L)==len(U),'Certificate domain')
    orders={(i,j) for i in range(n) for j in range(n) if U[i] is not None and U[i]<L[j]}
    need(cert['exact_checked'] is True and cert['forced_orders']==len(orders),'Certificate order metadata')
    if cert['status']=='feasible':
        beta=list(map(Q,cert['beta']));y=list(map(Q,cert['response']))
        need(len(beta)==len(X[0]) and len(y)==n,'Primal dimensions')
        z=[sum(a*b for a,b in zip(row,beta)) for row in X]
        need(z==list(map(Q,cert['latent'])),'Latent reconstruction')
        need(all(z[j]-z[i]>=1 for i,j in orders),'Exact unit margin')
        need(all(L[i]<=y[i] and (U[i] is None or y[i]<=U[i]) for i in range(n)),'Boundary response outside intervals')
        need(all(y[i]<=y[j] for i in range(n) for j in range(n) if z[i]<=z[j]),'Common nondecreasing response including ties')
        return {'status':'feasible','forced_orders':len(orders)}
    need(cert['status']=='infeasible','Unexpected certificate status')
    total=[Q(0)]*len(X[0]);weights=[];seen=set()
    for item in cert['dual']:
        i,j=item['lower_vertex'],item['higher_vertex'];weight=Q(item['weight'])
        need((i,j) in orders and (i,j) not in seen and weight>0,'Invalid positive dual support')
        seen.add((i,j));row=[b-a for a,b in zip(X[i],X[j])]
        need(row==item['difference'],'Stored feature difference mismatch')
        total=[a+weight*b for a,b in zip(total,row)];weights.append(weight)
    need(weights and sum(weights)==1 and all(x==0 for x in total),'Exact positive dual cancellation')
    return {'status':'infeasible','forced_orders':len(orders),'support_size':len(weights)}

def order_expiries(L,U,mask,metric):
    initial=Q(0) if metric=='additive' else Q(1)
    events={initial};lifetimes={}
    for i,u in enumerate(U):
        if u is None:continue
        for j,l in enumerate(L):
            if l<=u:continue
            degree=mask[i]+mask[j]
            if degree==0 or metric=='multiplicative_squared' and u==0:
                expiry=None
            elif metric=='additive':expiry=(l-u)/degree
            elif degree==1:expiry=l*l/(u*u)
            else:expiry=l/u
            lifetimes[i,j]=expiry
            if expiry is not None:events.add(expiry)
    return sorted(events),lifetimes

def expanded(L,U,mask,t,metric):
    if metric=='additive':
        return [max(Q(0),l-mask[i]*t) for i,l in enumerate(L)], [None if u is None else u+mask[i]*t for i,u in enumerate(U)]
    return [l*l/(t if mask[i] else 1) for i,l in enumerate(L)], [None if u is None else u*u*(t if mask[i] else 1) for i,u in enumerate(U)]

def verify_radius(record,X,L,U,mask,metric):
    need(record['metric']==metric,'Radius metric label')
    events,lifetimes=order_expiries(L,U,mask,metric)
    need(record['breakpoints']==list(map(str,events)),'Complete exact event inventory')
    tested={}
    for text,certificate in record['tested'].items():
        k=int(text);need(str(k)==text and 0<=k<len(events),'Invalid tested event index')
        lo,hi=expanded(L,U,mask,events[k],metric)
        tested[k]=certify(certificate,X,lo,hi)
    need(0 in tested,'Initial parameter not checked')
    if record['status']=='infinite':
        need(len(events)-1 in tested,'Permanent system not checked')
        lo,hi=expanded(L,U,mask,events[-1],metric)
        cert=record['permanent_dual'];certify(cert,X,lo,hi)
        need(cert['status']=='infeasible','Infinite result has no contradiction')
        need(all(lifetimes[d['lower_vertex'],d['higher_vertex']] is None for d in cert['dual']),'Nonpermanent infinity support')
        return {'status':'infinite','tested_certificates':len(tested),'breakpoints':len(events)}
    need(record['status']=='finite','Unknown radius status')
    t=Q(record['parameter_exact']);need(t in events,'Boundary is not a candidate event')
    at=events.index(t);need(at in tested and tested[at]['status']=='feasible','Boundary search evidence')
    need(record['parameter_meaning']==('epsilon' if metric=='additive' else 'gamma_squared'),'Parameter units')
    lo,hi=expanded(L,U,mask,t,metric)
    need(record['boundary_lower']==list(map(str,lo)) and record['boundary_upper']==[None if x is None else str(x) for x in hi],'Boundary interval reconstruction')
    need(lo[46]==hi[46]==1,'Native normalization anchor moved')
    upper=record['feasible_at_boundary'];need(upper['status']=='feasible','Missing upper witness');certify(upper,X,lo,hi)
    need(upper==record['tested'][str(at)],'Boundary certificate metadata disagrees')
    minimum=Q(0) if metric=='additive' else Q(1)
    support=[]
    if t==minimum:
        need(record['sharp_lower_certificate'] is None,'Initial feasible result has spurious lower domain')
    else:
        lower=record['sharp_lower_certificate'];need(lower and lower['status']=='infeasible','Missing sharp lower dual')
        need(at>0 and at-1 in tested,'Immediately preceding event absent')
        a,b=expanded(L,U,mask,events[at-1],metric);certify(lower,X,a,b)
        need(lower==record['tested'][str(at-1)],'Lower certificate event mismatch')
        for d in lower['dual']:
            i,j=d['lower_vertex'],d['higher_vertex'];expiry=lifetimes[i,j]
            need(expiry is None or expiry>=t,'A supported strict order expires below the claimed radius')
            support.append({'lower_vertex':i,'higher_vertex':j,'weight':d['weight'],
                            'expiry':None if expiry is None else str(expiry),'persists_for_all_smaller_parameters':True})
    meaning='original units' if metric=='additive' else 'square roots of returned rational responses'
    need(record['response_meaning']==meaning,'Response transformation interpretation')
    gamma=None
    if metric=='multiplicative_squared':
        gamma=record['gamma_enclosure'];a,b=Q(gamma['lower']),Q(gamma['upper'])
        need(a*a<=t<=b*b and a<=b and b-a<=Q(1,10**12),'Rational gamma enclosure')
        need(Q(gamma['decimal_lower_exact'])==a and Q(gamma['decimal_upper_exact'])==b,'Exact decimal gamma enclosure')
    return {'status':'finite','parameter_exact':str(t),'tested_certificates':len(tested),
            'breakpoints':len(events),'sharp_lower_support':support,'gamma_enclosure':gamma}

def run(primary,result,contract):
    need(result['anchor_vertex']==46 and result['source_sha256']==PRIMARY_SHA,'Anchor or primary attribution')
    rows=primary['vertices_data'];sequences=contract['sequences']
    need(len(rows)==len(sequences)==72 and len(set(sequences))==72,'Full72 domain')
    need([row['vertex'] for row in rows]==list(range(72)) and [row['sequence'] for row in rows]==sequences,'Identity order')
    need(sequences[46]==WT,'Native identity index')
    X,states=design(sequences);need([row['codon_states'] for row in rows]==states,'Sequence-to-codon reconstruction')
    values=[]
    for row in rows:
        need(len(row['scenario_ratios'])==27 and all(v is not None for v in row['scenario_ratios']),'Complete27 joint fields')
        v=list(map(Q,row['scenario_ratios']));need(min(v)>=0,'Negative relative count field');values.append(v)
    L=[min(v) for v in values];U=[max(v) for v in values]
    need(L==[Q(r['lower']) for r in rows] and U==[Q(r['upper']) for r in rows],'Primary empirical envelope')
    need(L[46]==U[46]==1,'Fixed reference interval')
    mask=[int(i!=46) for i in range(72)]
    need(set(result['full'])=={'additive','multiplicative_squared'},'Full metric coverage')
    full={metric:verify_radius(result['full'][metric],X,L,U,mask,metric) for metric in ('additive','multiplicative_squared')}
    deletions=result['all_single_vertex_deletions'];need(len(deletions)==72,'All72 row deletions required')
    row_summaries=[]
    for omitted,row in enumerate(deletions):
        ids=[v for v in range(72) if v!=omitted]
        need(row['omitted_vertex']==omitted and row['retained_vertex_ids']==ids,'Deleted-row identity or coverage')
        c=certify(row['certificate'],[X[i] for i in ids],[L[i] for i in ids],[U[i] for i in ids])
        row_summaries.append({'omitted_vertex':omitted,'kind':'fixed_normalization_identity' if omitted==46 else 'measured_genotype_observation',**c})
    prep=result['all_single_preparation_deletions'];need(len(prep)==12,'All12 preparation deletions required')
    prep_summaries=[]
    for index,row in enumerate(prep):
        kind='output' if index<9 else 'input';omitted=index if index<9 else index-9
        columns=[3*b+t for b in range(9) for t in range(3) if (b if kind=='output' else t)!=omitted]
        need(row['kind']==kind and row['omitted_zero_based']==omitted and row['retained_scenarios']==columns,'Preparation selection or ordering')
        lo=[min(v[k] for k in columns) for v in values];hi=[max(v[k] for k in columns) for v in values]
        need(all(a<=b<=c<=d for a,b,c,d in zip(L,lo,hi,U)),'Deleted-scenario envelope nesting')
        need(set(row['results'])=={'additive','multiplicative_squared'},'Preparation metric coverage')
        checked={metric:verify_radius(row['results'][metric],X,lo,hi,mask,metric) for metric in ('additive','multiplicative_squared')}
        for metric,radius in checked.items():
            if full[metric]['status']=='infinite':need(radius['status']=='infinite','Narrower set repaired infinite exclusion')
            elif radius['status']=='finite':need(Q(radius['parameter_exact'])>=Q(full[metric]['parameter_exact']),'Nesting reduced radius')
        prep_summaries.append({'kind':kind,'omitted_zero_based':omitted,'retained_scenario_count':len(columns),'radii':checked})
    all_radii=list(full.values())+[v for row in prep_summaries for v in row['radii'].values()]
    feasible=[r['omitted_vertex'] for r in row_summaries if r['status']=='feasible']
    return {'status':'PASS','radii_checked':len(all_radii),'row_deletion_certificates_checked':72,
            'preparation_omission_sets_checked':12,'searched_event_certificates_checked':sum(r['tested_certificates'] for r in all_radii),
            'complete_directed_event_inventories_checked':26,'full':full,'row_deletions':row_summaries,
            'preparation_deletions':prep_summaries,'feasible_row_omissions':feasible,
            'feasible_measured_observation_omissions':[v for v in feasible if v!=46],
            'reference_identity_omission_feasible':46 in feasible,
            'maximum_exact_discrepancy':'0','reoptimization_performed':False,
            'interpretation':'Reference46 is a normalization identity fixed at1. Deleting its constraint is structurally different from deleting a measured genotype observation.',
            'scope':'Independent exact certificates and reconstruction of frozen exploratory outputs, not biological calibration or new validation.'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repo','--root',dest='repo',type=Path,required=True)
    parser.add_argument('--result',type=Path);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    result_path=args.result if args.result is not None else args.repo/'results/cycle7/robustness.json.gz'
    need(not args.out.exists(),'Refusing existing output')
    start=time.perf_counter();started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    frozen_path=args.repo/'analysis/cycle7/freeze.json';frozen_raw=read_pinned(frozen_path,FREEZE_SHA);frozen=json.loads(frozen_raw)
    input_hashes={str(frozen_path):FREEZE_SHA}
    for name,digest in frozen['files'].items():
        path=frozen_path.parent/name;read_pinned(path,digest);input_hashes[str(path)]=digest
    core=args.repo/'analysis/cycle6/fas_pilot/ordinal_geometry.py';read_pinned(core,frozen['ordinal_core_sha256']);input_hashes[str(core)]=frozen['ordinal_core_sha256']
    saved_result=result_path.read_bytes()
    result_raw=gzip.decompress(saved_result) if result_path.suffix=='.gz' else saved_result
    need(sha(result_raw)==RESULT_SHA,'Decoded successor result hash mismatch')
    input_hashes[str(result_path)]=sha(saved_result)
    primary_path=args.repo/'results/cycle6/full_result_v1.json.gz';compressed=read_pinned(primary_path,GZIP_SHA)
    input_hashes[str(primary_path)]=GZIP_SHA;raw=gzip.decompress(compressed);need(sha(raw)==PRIMARY_SHA,'Decompressed primary hash')
    contract_path=args.repo/'analysis/cycle6/fas_pilot/source_contract.json';contract_raw=read_pinned(contract_path,CONTRACT_SHA)
    input_hashes[str(contract_path)]=CONTRACT_SHA
    audit=run(json.loads(raw),json.loads(result_raw),json.loads(contract_raw))
    for path,digest in input_hashes.items():read_pinned(Path(path),digest)
    audit.update({'input_hashes_before_and_after':input_hashes,'decompressed_primary_sha256':PRIMARY_SHA,
                  'verifier_sha256':sha(Path(__file__).read_bytes()),'started_at':started,
                  'runtime_seconds':time.perf_counter()-start,'python':sys.version,'platform':platform.platform()})
    args.out.mkdir(parents=True)
    (args.out/'receipt.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:audit[key] for key in ('status','radii_checked','row_deletion_certificates_checked','searched_event_certificates_checked',
                                               'feasible_row_omissions','maximum_exact_discrepancy','runtime_seconds')}))

if __name__=='__main__':main()
