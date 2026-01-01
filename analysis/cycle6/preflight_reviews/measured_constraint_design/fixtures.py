"""Small exact structural counterexamples. No biological data or author imports."""
if not __debug__:raise RuntimeError('Optimization refused')
from itertools import product
from collections import deque,Counter
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def graph(domains):
 vertices=list(product(*[range(n) for n in domains]))
 adjacency={u:[v for v in vertices if sum(a!=b for a,b in zip(u,v))==1] for u in vertices}
 return vertices,adjacency
def distance(adj,allowed,s,t):
 if s not in allowed or t not in allowed:return None
 q=deque([(s,0)]);seen={s}
 while q:
  u,d=q.popleft()
  if u==t:return d
  for v in adj[u]:
   if v in allowed and v not in seen:seen.add(v);q.append((v,d+1))
 return None
# Complete standard-code three-Lys synonymous geometry, AAA/AAG at each codon.
v,g=graph([2,2,2]);s=(0,0,0);t=(1,1,0)
F={(1,0,0),(0,1,0)};D={x for x in v if sum(x)<=2}
require(distance(g,set(v),s,t)==2,'Baseline distance')
require(distance(g,D-F,s,t) is None,'Bounded-domain separator')
require(distance(g,set(v)-F,s,t)==4,'Unmeasured admissible bypass not preserved')
unknown=sorted(set(v)-F-{s,t});completion_distances=Counter()
for bits in product([False,True],repeat=len(unknown)):
 allowed={s,t}|{x for x,b in zip(unknown,bits) if b}
 d=distance(g,allowed,s,t)
 require(d is None or d>=2*(2+1)-2,'Escape lower bound violated')
 completion_distances[str(d)]+=1
require(dict(completion_distances)=={'None':13,'4':3},'Exhaustive completion count')
# A third allowed allele can bypass a sampled Q2 without a third position.
v2,g2=graph([3,2]);s2=(0,0);t2=(1,1);F2={(1,0),(0,1)}
require(distance(g2,set(v2)-F2,s2,t2)==3,'Third-allele bypass')
# One unresolved cut vertex destroys the exclusion. Weak >= admits ties.
require(distance(g,set(v)-{(1,0,0)},s,t)==2,'Unknown is not failure')
tau=Q(1,2)
intervals={s:(Q(3,4),Q(1)),t:(Q(3,4),Q(1)),(1,0,0):(Q(0),tau),(0,1,0):(Q(0),Q(1,4))}
possible={x for x in v if x not in intervals or intervals[x][1]>=tau}
require(distance(g,possible&D,s,t)==2,'Tie incorrectly excluded')
# Joint-only obstruction can be exactly additive in each trait.
vs,gs=graph([2,2]);a=(0,0);b=(1,1)
f1={x:Q(x[0]-x[1]) for x in vs};f2={x:-f1[x] for x in vs}
require(distance(gs,{x for x in vs if f1[x]>=0},a,b)==2,'Trait1 path')
require(distance(gs,{x for x in vs if f2[x]>=0},a,b)==2,'Trait2 path')
require(distance(gs,{x for x in vs if f1[x]>=0 and f2[x]>=0},a,b) is None,'Joint additive obstruction')
# Scalar reciprocal sign pattern has strictly positive interaction contrast.
scalar={a:Q(3,4),b:Q(3,4),(1,0):Q(1,4),(0,1):Q(1,4)}
contrast=scalar[a]+scalar[b]-scalar[(1,0)]-scalar[(0,1)]
require(contrast==1,'Exact scalar contrast')
# Bonferroni can be sharp under dependent marginal failures: four disjoint
# failure events, each with probability1/20, leave simultaneous coverage4/5.
worlds=range(20);failure_events=[{i} for i in range(4)]
union=set().union(*failure_events)
require(Q(len(union),20)==sum((Q(len(e),20) for e in failure_events),Q(0))==Q(1,5),'Union bound fixture')
# Parent-supplied FAS structural type: arbitrary codon unary scores on P3
# can have a gate valley. This is not an independent sequence reconstruction.
vp=list(product(range(3),range(3),range(2),range(2),range(2)))
gp={u:[w for w in vp if sum(a!=b for a,b in zip(u,w))==1 and sum(abs(a-b) for a,b in zip(u,w))==1] for u in vp}
sp=(0,0,0,0,0);tp=(2,0,0,0,0)
require(len(vp)==72 and sum(map(len,gp.values()))//2==204,'Product structural count')
unary={u:Q(0) if u[0]==1 else Q(1) for u in vp}
require(distance(gp,set(vp),sp,tp)==2 and distance(gp,{u for u in vp if unary[u]>=Q(1,2)},sp,tp) is None,'Codon-unary gate valley')
require(1+sum(n-1 for n in [3,3,2,2,2])==8,'Codon-unary model dimension')
result={'status':'EXACT_SYNTHETIC_FIXTURES_PASS','Q3_bounded_domain_radius':2,'Q3_baseline_distance':2,'Q3_certified_detour_minimum_if_route_exists':4,'Q3_unmeasured_completion_counts':dict(completion_distances),'Q3_global_connection_unresolved':True,'K3xK2_third_allele_bypass_length':3,'unknown_and_threshold_tie_kept_possible':True,'joint_only_additive_counterexample_verified':True,'scalar_interaction_contrast':'1','four95percent_marginals_can_have_simultaneous_coverage':'4/5','P3xP3xQ3_vertices':72,'P3xP3xQ3_edges':204,'codon_unary_scalar_disconnection_without_intercodon_interaction':True,'codon_unary_dimension_including_intercept':8,'scope':'Finite invented structural examples, not proof of statistical calibration, independent FAS sequence reconstruction or biological results','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out=Path(__file__).with_name('fixtures_v2.json')
if out.exists():raise RuntimeError('No overwrite')
out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
