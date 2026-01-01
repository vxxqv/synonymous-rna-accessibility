"""Post-freeze schema adapter. Imports only the independently frozen verifier."""
from pathlib import Path
import collections, gzip, hashlib, itertools, json, math, sys, time
import numpy as np
sys.dont_write_bytecode=True
import verifier as V
H=Path(__file__).resolve().parent; R=H.parents[2]; P=R/'work/gill_multiclass_cycle3'
def read(p):
 b=Path(p).read_bytes();return json.loads(gzip.decompress(b) if str(p).endswith('.gz') else b)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
pins={'validation_v1/receipt.json':'5f7535bb58622cb656ee0cb7ed4109f8daef324ac470eb948ba40ffd9ab6757c','validation_v1/class_metrics.json.gz':'f1c4c1fac09362edda68aa20b51c386784905837ca3eb9f89f256c7148e95491','validation_v1/heldout_counts.json.gz':'bcbf8989413b6fa5c9c1a3c945c75693f5b4a0e194a857754a41cc61848f09c9','validation_v1/validation.json':'cdbdd83d6861c54d1cd5094990cbeb306803349c63bdd0519caf391676c819f6','training_v1/receipt.json':'c6fb3309a93635efab38f1469516617c3f1b4f857c88b4fbfff91be21115bdef'}
for p,s in pins.items():assert sha(P/p)==s
assert sha(H/'verifier.py')=='4331bb0ed7349bfc896f59bdbdfb98b931818c6671b9fa908c8a261a5cdfa588'
assert sha(H/'evaluation_run_v1/calculation/independent_results.json.gz')=='b0a019cf1754598b672fba59af497ec9aefd1effabfa1515bb4cd15e695a2a82'
I=read(H/'evaluation_run_v1/calculation/independent_results.json.gz');A=read(P/'validation_v1/class_metrics.json.gz');S=read(P/'validation_v1/validation.json');T=read(P/'training_v1/predictions.json.gz');C=read(P/'training_v1/counts.json.gz');D=read(P/'validation_v1/heldout_counts.json.gz')
inventory=read(V.BASE/'inventory.json.gz'); meta={c['class_id']:c for c in inventory['classes']}
fail=[]; counts=collections.Counter();maxerr=collections.defaultdict(float); examples={}; numerical_notes=[]
def exact(cat,x,y,loc):
 counts[cat]+=1
 if x!=y:fail.append({'category':cat,'location':loc,'independent':x,'author':y})
def close(cat,x,y,loc,tol=1e-9):
 counts[cat]+=1
 if x is None or y is None:
  if x!=y:fail.append({'category':cat,'location':loc,'independent':x,'author':y,'reason':'undefined_status'})
  return
 x=float(x);y=float(y);e=abs(x-y)
 if e>maxerr[cat]:maxerr[cat]=e;examples[cat]=loc
 if not math.isfinite(e) or e>tol:fail.append({'category':cat,'location':loc,'independent':x,'author':y,'error':e,'tolerance':float(tol)})
def array(cat,x,y,loc,tol):
 x=np.asarray(x);y=np.asarray(y);exact(cat+'_shape',list(x.shape),list(y.shape),loc)
 if x.shape!=y.shape:return
 for idx in np.ndindex(x.shape):close(cat,x[idx],y[idx],loc+':'+str(idx),tol)
pair3=list(itertools.combinations(range(3),2));start=time.monotonic();headlines=[]
exact('scenario_order',S['scenario_order'],[list(x) for x in V.SCENARIOS],'global')
exact('decomposition_order',S['decomposition_pair_order'],[list(x) for x in V.PAIR45],'global')
for cp in ['0.25','0.5','1.0']:
 for part in ['development','validation']:
  own=I[cp]['classes'][part];auth={r['class_id']:r for r in A[cp][part]};agg=I[cp]['aggregate'][part];summary=S['all_summaries'][cp][part]
  exact('membership',sorted(r['class_id'] for r in own),sorted(auth),cp+part)
  rankrows={'RNA':[],'DNA':[]};authorbarriers=collections.defaultdict(list)
  for row in own:
   cid=row['class_id'];a=auth[cid];n=row['n'];loc=f'{cp}/{part}/{cid}';m=meta[cid];keys,X,adj=V.setup_class(m)
   exact('vertices',n,a['vertices'],loc);exact('stratum',row['stratum'],a['stratum'],loc);exact('clique_flag',row['all_factor_cliques'],a['clique_product'],loc)
   exact('rank',row['stats']['projection_rank'],a['codon_rank'],loc);exact('residual_dimension',row['stats']['residual_dimension'],a['residual_dimensions'],loc)
   exact('training_keys',keys,T[cid]['keys'],loc);exact('training_partition',part,T[cid]['partition'],loc)
   delta=256*np.finfo(float).eps*row['raw_log_input_norm']+1e-10*row['centered_log_input_norm']
   dot_tol=delta*(4*row['raw_log_input_norm']+2*delta)
   for ok,ak in [('residual_cosine','cosine'),('residual_mse','residual_prediction_mse'),('zero_residual_mse','zero_residual_mse'),('whole_cosine','whole_centered_cosine'),('whole_mse','whole_prediction_mse'),('codon_main_prediction_mse','codon_prediction_mse'),('constant_prediction_mse','constant_prediction_mse')]:
    close('class_'+ak,row['stats'][ok],a['stats'][ak],loc,1e-9 if 'cosine' in ok else dot_tol/n)
   for ok,ak in [('residual_TT','train_energy'),('residual_TV','cross'),('residual_VV','validation_energy')]:close('class_'+ak,row['stats'][ok]/n,a['stats'][ak],loc,dot_tol/n)
   array('primary_terms',row['primary_terms']+[row['stats']['residual_TV']],a['primary_decomposition_terms'],loc,dot_tol)
   array('decomposition_terms',[d['terms']+[d['total']] for d in row['decomposition45']],a['decomposition_terms'],loc,dot_tol)
   for kind in ['RNA','DNA']:array('gram',np.asarray(row[kind+'_residual_gram'])/n,a[kind+'_residual_gram'],loc+kind,dot_tol/n)
   for k,(i,j) in enumerate(pair3):
    q=row['pair_stats'][str((i,j))]
    for kind in ['RNA','DNA']:close('residual_pair_cosine',q[kind+'_residual_cosine'],a[kind+'_residual_cosines'][k],loc+kind+str(k))
    close('RNA_rank',q['RNA_midrank_correlation'],a['RNA_rank_correlations'][k],loc+str(k))
   dna=np.array([[C[k][x] for x in ['E','F','G']] for k in keys]);dnrank=[V.cosine(V.midranks(dna[:,i]),V.midranks(dna[:,j])) for i,j in pair3]
   for k,x in enumerate(dnrank):close('DNA_rank',x,a['DNA_rank_correlations'][k],loc+str(k))
   rankrows['DNA'].append(dnrank);rankrows['RNA'].append([row['pair_stats'][str(p)]['RNA_midrank_correlation'] for p in pair3])
   tf=T[cid]['fields'][cp]
   for f in ['T','PT','QT']:
    av=np.asarray(tf[f]);iv=np.asarray(row['independent_training_fields'][f]);err=float(np.linalg.norm(av-iv));close('training_vector_norm',err,0,loc+f,delta)
   raw=np.array([[C[k]['B'],C[k]['C'],D[k]['D'],C[k]['E'],C[k]['F'],C[k]['G']] for k in keys],dtype=float)
   logs=np.log2(raw+float(cp));proj,res=V.project(X,logs)
   for kind,cols in [('RNA',[0,1]),('DNA',[3,4,5])]:
    exact('training_raw_counts',raw[:,cols].T.astype(int).tolist(),tf[kind+'_raw_counts'],loc+kind)
    array('training_log_fields',logs[:,cols].T,tf[kind+'_logs'],loc+kind,1e-12)
    for j,col in enumerate(cols):close('training_residual_vector_norm',float(np.linalg.norm(res[:,col]-np.asarray(tf[kind+'_residuals'][j]))),0,loc+kind+str(j),delta)
   scenarios=np.array([logs[:,i]-logs[:,3+j] for i,j in V.SCENARIOS]);ps=np.array([proj[:,i]-proj[:,3+j] for i,j in V.SCENARIOS])
   fields=np.vstack([scenarios,ps,np.array(row['independent_training_fields']['T'])[None,:],np.array(row['independent_training_fields']['PT'])[None,:]])
   # Reconstruct every all-pair matrix so every supplied maximum witness can be checked, independently of the author method.
   em=np.minimum(fields[:,:,None],fields[:,None,:]);cap=np.where((adj|np.eye(n,dtype=bool))[None,:,:],em,-np.inf)
   for k in range(n):cap=np.maximum(cap,np.minimum(cap[:,:,k,None],cap[:,None,k,:]))
   barrier=em-cap
   objs=a['scenario_barriers']+a['scenario_codon_barriers']+[tf['training_barrier'],tf['training_codon_barrier']]
   vals=row['scenario_barriers']+row['scenario_codon_barriers']+[row['training_barrier'],row['training_codon_barrier']]
   for j,(obj,value,F) in enumerate(zip(objs,vals,fields)):
    bl=loc+':field'+str(j);tol=1e-9*max(float(np.ptp(F)),np.finfo(float).eps*float(np.max(np.abs(F))),np.finfo(float).tiny)
    close('barrier_maximum',value,obj['maximum'],bl,tol);exact('endpoint_pair_count',n*(n-1)//2,obj['endpoint_pairs_covered'],bl)
    precision=128*np.finfo(float).eps*n*max(1,float(np.max(np.abs(F))))
    close('barrier_precision',precision,obj['precision_tolerance'],bl,1e-25+1e-12*precision)
    exact('barrier_positive_flag',obj['maximum']>obj['precision_tolerance'],obj['numerically_resolved_positive'],bl)
    exact('barrier_tiny_flag',0<obj['maximum']<=obj['precision_tolerance'],obj['unresolved_tiny_positive'],bl)
    pair=obj['arg_pair']
    if pair is None:exact('barrier_empty_witness',obj['maximum'],0,bl)
    else:
     exact('barrier_witness_valid',True,len(pair)==2 and all(isinstance(x,int) and 0<=x<n for x in pair) and pair[0]!=pair[1],bl)
     close('barrier_witness',float(barrier[j,pair[0],pair[1]]),obj['maximum'],bl,tol)
   for field,j in [('heldout_barrier',8),('heldout_codon_barrier',17)]:exact('heldout_barrier_alias',a[field],objs[j],loc+field)
   flag='clique' if row['all_factor_cliques'] else 'nonclique'
   authorbarriers[flag].append(a)
  sloc=f'{cp}/{part}'
  exact('aggregate_class_count',len(own),summary['class_count'],sloc)
  exact('target_value',0.5,summary['target'],sloc)
  exact('failure_status','target_passed' if agg['operational_R_ge_0_5'] else 'target_failed',summary['status'],sloc)
  exact('zero_energy_counts',{'training':sum(r['stats']['residual_TT']==0 for r in own),'validation':sum(r['stats']['residual_VV']==0 for r in own)},summary['zero_energy_classes'],sloc)
  exact('clique_class_count',sum(r['all_factor_cliques'] for r in own),summary['clique_classes'],sloc)
  exact('nonclique_class_count',sum(not r['all_factor_cliques'] for r in own),summary['nonclique_classes'],sloc)
  for ok,ak in [('residual_TT','train_energy'),('residual_TV','cross'),('residual_VV','validation_energy')]:close('aggregate_energy',agg['terms'][ok],summary[ak],sloc+ak,1e-9)
  close('aggregate_cosine',agg['pooled_residual_cosine'],summary['cosine'],sloc)
  exact('aggregate_target',agg['operational_R_ge_0_5'],summary['cosine']>=summary['target'],sloc)
  for ok,ak in [('residual_mse','residual_prediction_mse'),('zero_residual_mse','zero_residual_mse'),('whole_mse','whole_prediction_mse'),('codon_main_prediction_mse','codon_prediction_mse'),('constant_prediction_mse','constant_prediction_mse')]:close('aggregate_loss',agg['sum_class_mean_losses'][ok]/len(own),summary['losses'][ak],sloc+ak)
  array('aggregate_primary',agg['primary_decomposition_terms']+[agg['terms']['residual_TV']],summary['primary_decomposition_aggregate_terms'],sloc,1e-9)
  array('aggregate_decomposition',[x['terms']+[x['total']] for x in agg['decomposition45']],summary['decomposition_aggregate_terms'],sloc,1e-9)
  array('aggregate_scenario_cosine',[x['cosine'] for x in agg['decomposition45']],summary['scenario_residual_pooled_cosines'],sloc,1e-9)
  for kind in ['RNA','DNA']:
   vals=[agg['RNA_DNA_pooled_pair_cosines'][kind+'_residual_gram'][str(p)] for p in pair3]
   array('aggregate_library_cosine',vals,summary[kind+'_residual_pooled_cosines'],sloc+kind,1e-9)
   for j in range(3):
    r=[x[j] for x in rankrows[kind] if x[j] is not None];su=summary['rank_agreement_summary'][kind+'_rank_correlations'][j]
    exact('rank_summary_defined',len(r),su['defined'],sloc+kind+str(j));exact('rank_summary_undefined',len(own)-len(r),su['undefined'],sloc+kind+str(j))
    close('rank_summary_mean',math.fsum(r)/len(r) if r else None,su['mean_defined'],sloc+kind+str(j));array('rank_summary_range',[min(r),max(r)] if r else [],su['range_defined'],sloc+kind+str(j),1e-9)
  alo={x['stratum']:x for x in summary['leave_one_stratum_out']};ilo={x['stratum']:x for x in agg['leave_one_stratum_out']};exact('influence_membership',sorted(ilo),sorted(alo),sloc)
  for st,x in ilo.items():
   close('influence_cosine',x['pooled_residual_cosine_without'],alo[st]['remaining_cosine'],sloc+st)
   rows=[r for r in own if r['stratum']==st]
   for ik,ak in [('residual_TT','train_energy_fraction'),('residual_VV','validation_energy_fraction')]:close('influence_energy',math.fsum(r['stats'][ik]/r['n'] for r in rows)/agg['terms'][ik],alo[st][ak],sloc+st+ik)
  array('influence_range',agg['leave_one_stratum_out_range'],summary['leave_one_stratum_out_range'],sloc,1e-9)
  for ik,ak in [('residual_TT','largest_stratum_training_energy_fraction'),('residual_VV','largest_stratum_validation_energy_fraction')]:close('influence_maximum',agg['maximum_stratum_energy_fraction'][ik],summary[ak],sloc+ik)
  for flag,rows in authorbarriers.items():
   sb=summary['barriers_by_architecture'][flag];exact('barrier_summary_classes',len(rows),sb['classes'],sloc+flag)
   for key in ['scenario_barriers','scenario_codon_barriers']:
    for j in range(9):
     values=[r[key][j] for r in rows];sbj=sb[key][j]
     close('barrier_summary_max',max(v['maximum'] for v in values),sbj['maximum'],sloc+flag+key+str(j))
     close('barrier_summary_mean',math.fsum(v['maximum'] for v in values)/len(values),sbj['mean'],sloc+flag+key+str(j))
     for out,fn in [('numerically_resolved_positive_classes',lambda v:v['numerically_resolved_positive']),('raw_positive_classes',lambda v:v['maximum']>0),('unresolved_tiny_positive_classes',lambda v:v['unresolved_tiny_positive'])]:exact('barrier_summary_counts',sum(fn(v) for v in values),sbj[out],sloc+flag+key+str(j)+out)
  headlines.append({'pseudocount':cp,'partition':part,'classes':len(own),'vertices':agg['vertices'],'strata':agg['strata'],'residual_cosine':agg['pooled_residual_cosine'],'target_pass':agg['operational_R_ge_0_5'],'residual_prediction_gain_vs_zero':agg['residual_loss_gain_vs_zero'],'leave_one_stratum_range':agg['leave_one_stratum_out_range']})
exact('primary_summary_alias',S['all_summaries']['0.5']['validation'],S['primary_result'],'global')
exact('retained_failed_pilot',True,S['retained_failed_reference_pilot'],'global')
out={'status':'pass' if not fail else 'discrepancies','checks':dict(counts),'maximum_absolute_deviations':dict(maxerr),'maximum_deviation_locations':examples,'discrepancies':fail,'headlines':headlines,'input_hashes':pins,'source_sha256':sha(__file__),'frozen_verifier_sha256':sha(H/'verifier.py'),'runtime_seconds':time.monotonic()-start,'scope':'All2176classes,3pseudocounts,45normalizationpairs,130560fieldgraphs and supplied maximum witnesses. No workbook access, no author implementation import.'}
(H/'comparison_v3.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':out['status'],'discrepancies':len(fail),'checks':sum(counts.values()),'runtime':out['runtime_seconds'],'first_discrepancies':fail[:8]},indent=2))
