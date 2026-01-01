"""Read-only finite corroboration using threshold connectivity, not heap paths."""
import collections,hashlib,importlib.util,itertools,json,math,pathlib,random
REPO=pathlib.Path(__file__).resolve().parents[2]
if not __debug__:
 raise RuntimeError('Exact verification requires enabled assertions')
SOURCE=REPO/'src/barrier_identification.py'
spec=importlib.util.spec_from_file_location('audited_barrier_source',SOURCE)
api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
counts=collections.Counter()
def connected(adj,retained,a,b):
 if a not in retained or b not in retained:return False
 seen={a};todo=[a]
 while todo:
  for v in adj[todo.pop()]:
   if v in retained and v not in seen:seen.add(v);todo.append(v)
 return b in seen
def bottleneck(adj,f,a,b):
 for t in sorted(set(f),reverse=True):
  if connected(adj,{v for v,x in enumerate(f) if x>=t},a,b):return t
 raise ValueError('disconnected')
def reference(adj,f,a,b):return min(f[a],f[b])-bottleneck(adj,f,a,b)
def all_graphs(n):
 pairs=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  adj=[[]for _ in range(n)]
  for i,(a,b) in enumerate(pairs):
   if mask>>i&1:adj[a].append(b);adj[b].append(a)
  yield adj
def product_graph(g,h):
 n=len(h)
 return [[y*n+z for y in g[x]]+[x*n+w for w in h[z]] for x in range(len(g)) for z in range(n)]
rng=random.Random(15112026)
for n in range(1,5):
 for adj in all_graphs(n):
  counts['all_labelled_graphs_n1_to4']+=1
  for a in range(n):assert api.barrier(adj,[0]*n,a,a)==0;counts['coincident_endpoint_checks']+=1
  for a,b in itertools.combinations(range(n),2):
   if not connected(adj,set(range(n)),a,b):
    for call in (lambda:api.barrier(adj,[0]*n,a,b),lambda:api.sharp_barrier_bounds(adj,[0]*n,[1]*n,a,b)):
     try:call()
     except ValueError:counts['disconnected_endpoint_rejections']+=1
     else:raise AssertionError('disconnected baseline accepted')
    continue
   for f in itertools.product((-1,0,1),repeat=n):
    expected=reference(adj,f,a,b);assert api.barrier(adj,f,a,b)==expected
    counts['field_comparisons_threshold_oracle']+=1
    phi=tuple(x*x*x+2*x for x in f)
    bb=bottleneck(adj,f,a,b);mm=min(f[a],f[b])
    assert api.barrier(adj,phi,a,b)==(mm**3+2*mm)-(bb**3+2*bb)
    assert expected>0 or api.barrier(adj,tuple(max(0,x)for x in f),a,b)==0
    counts['transform_checks']+=2
   for _ in range(12):
    box=[tuple(sorted(rng.choices((-2,-1,0,2),k=2)))for _ in range(n)]
    lower,upper=zip(*box)
    vals=[reference(adj,f,a,b)for f in itertools.product(*(set(p)for p in box))]
    assert api.sharp_barrier_bounds(adj,lower,upper,a,b)==(min(vals),max(vals))
    counts['sampled_box_corner_ranges']+=1
    for t in (-1,0,2):
     if min(lower[a],lower[b])<t:expected='endpoints_not_certified'
     else:
      decisions=[connected(adj,{v for v,x in enumerate(f)if x>=t},a,b)for f in itertools.product(*(set(p)for p in box))]
      expected='connected'if all(decisions)else'separated'if not any(decisions)else'unresolved'
     assert api.threshold_status(adj,lower,upper,a,b,t)==expected
     counts['threshold_corner_decisions']+=1
g_fixtures=[[[1],[0,2],[1]],[[1,2],[0,2],[0,1]],[[1,3],[0,2],[1,3],[0,2]]]
h=[[1],[0,2],[1]]
for g in g_fixtures:
 gh=product_graph(g,h)
 for hv in itertools.product((-2,0,3),repeat=3):
  for z in range(3):
   hmax=max(hv)
   has_ascent=any(hv[w]==hmax and connected(h,{v for v,x in enumerate(hv)if x>=hv[z]},z,w)for w in range(3))
   if not has_ascent:continue
   for gv in itertools.product((-2,0,3),repeat=len(g)):
    f=[x+y for x in gv for y in hv]
    for a,b in itertools.combinations(range(len(g)),2):
     local=reference(g,gv,a,b)
     expected=max(0,local-(hmax-hv[z]))
     actual=reference(gh,f,3*a+z,3*b+z)
     assert actual==expected
     counts['arbitrary_background_ascent_product_cases']+=1
path=[[1],[0,2],[1]]
mixed_range={'lower':[0,0,0],'upper':[1,1,1],'barrier_at_all_lower':reference(path,[0,0,0],0,2),'barrier_at_all_upper':reference(path,[1,1,1],0,2),'sharp_range':list(api.sharp_barrier_bounds(path,[0,0,0],[1,1,1],0,2)),'correlated_equal_score_completion_range':[0,0]}
counter_g=[2,0,2];counter_h=[0,-10,3];prod=product_graph(path,path)
counteractual=reference(prod,[x+y for x in counter_g for y in counter_h],0,6)
assert counteractual==2
longpath=[[1],[0,2],[1,3],[2,4],[3]];gv=[2,20,0,20,2];hv=[0,-1,3]
notnecessary=reference(product_graph(longpath,path),[x+y for x in gv for y in hv],0,12)
assert notnecessary==0 and max(0,reference(longpath,gv,0,4)-3)==0
assert reference(path,[2,-1,2],0,2)-reference(path,[1,0,1],0,2)==2
nonfinite={}
for key,call in [('helper_nan_local',lambda:api.additive_background_barrier(float('nan'),0,1)),('threshold_nan',lambda:api.threshold_status(path,[0,0,0],[1,1,1],0,2,float('nan')))]:
 try:nonfinite[key]={'accepted':call()}
 except ValueError as exc:nonfinite[key]={'rejected_with_ValueError':str(exc)}
result={'status':'passed_for_stated_finite_inputs','scope':'Independent threshold-connectivity implementation; finite corroboration, not universal proof; no biological data','counts':dict(counts),'mixed_endpoint_internal_box_example':mixed_range,'missing_background_assumption_counterexample':{'g':counter_g,'h':counter_h,'actual_B':counteractual,'formula_without_hypothesis':0},'background_assumption_not_instancewise_necessary':{'g':gv,'h':hv,'actual_B':notnecessary,'formula':0},'perturbation_constant_2_attained':True,'nonfinite_out_of_domain_validation_observations':nonfinite,'sha256':{'implementation':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'audit_script':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}}
out=REPO/'results/programme/identification_v1/independent_replay.json';out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
