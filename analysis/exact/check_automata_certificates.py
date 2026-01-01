"""Standalone integer certificate checker; imports no research implementation."""
from pathlib import Path
from collections import defaultdict
import ast,json,hashlib
import networkx as nx
base=Path(__file__).resolve().parent
r=json.loads((base/"exact_peptide_automata.json").read_text())
out={"feasibility":{},"threshold":{}}
for d,z in r.items():
 M={int(s):{int(t):v for t,v in row.items()} for s,row in z["transition_matrix"].items()}
 cert=z.get("exact_pf_certificate")
 summary={"states":len(M),"max_row_sum":max(sum(row.values()) for row in M.values())}
 if cert:
  x={int(s):int(v) for s,v in cert["positive_integer_vector"].items()}
  assert set(x)==set(M) and min(x.values())>0
  lo,hi,den=[int(cert[k]) for k in ("lower_numerator","upper_numerator","denominator")]
  assert all(lo*x[s] <= den*sum(v*x[t] for t,v in row.items()) <= hi*x[s] for s,row in M.items())
  summary["rational_PF_enclosure_verified"]=[str(lo)+"/"+str(den),str(hi)+"/"+str(den)]
 def reachable(graph,start):
  seen={start};todo=[start]
  for s in todo:
   for t in graph[s]:
    if t not in seen:seen.add(t);todo.append(t)
  return seen
 s=next(iter(M));reverse={s:set() for s in M}
 for a,row in M.items():
  for b in row:reverse[b].add(a)
 assert len(reachable(M,s))==len(M)==len(reachable(reverse,s))
 assert any(s in row for s,row in M.items())
 summary["one_primitive_SCC_verified"]=True
 out["feasibility"][d]=summary
 current={int(s):int(c) for s,c in z['initial'].items()}
 for term in z['exact_counts']:
  assert sum(current.values())==int(term['both_nonempty'])
  following=defaultdict(int)
  for s,c in current.items():
   for t,v in M[s].items():following[t]+=c*v
  current=following
 summary['100_exact_count_terms_verified']=True
q=json.loads((base/"exact_bridge_threshold_automata.json").read_text())
for d,z in q.items():
 M={ast.literal_eval(s):{ast.literal_eval(t):v for t,v in row.items()} for s,row in z["matrix"].items()}
 current={ast.literal_eval(s):c for s,c in z["initial"].items()}
 vals=[]
 for m in range(1,101):
  n0=sum(c for s,c in current.items() if s[0])
  n1=sum(c for s,c in current.items() if not s[0])
  assert str(n0)==z["exact_counts"][m-1]["delta0_peptide_pairs"]
  assert str(n1)==z["exact_counts"][m-1]["delta1_peptide_pairs"]
  if d=="1":
   assert n0==24*9**(m-1)
   assert n1==8*(3*m-1)*9**(m-1)
  nxt=defaultdict(int)
  for s,c in current.items():
   for t,v in M[s].items():nxt[t]+=c*v
  current=nxt
  vals.append((n0,n1))
 checked=0
 graph=nx.DiGraph();graph.add_nodes_from(M)
 graph.add_edges_from((s,t) for s,row in M.items() for t in row)
 reachable_states=set()
 for s in z['initial']:
  state=ast.literal_eval(s);reachable_states.add(state);reachable_states.update(nx.descendants(graph,state))
 assert reachable_states==set(M)
 cyclic_components={frozenset(c) for c in nx.strongly_connected_components(graph) if len(c)>1 or any(graph.has_edge(s,s) for s in c)}
 certified_components=set()
 for c in z["nontrivial_SCCs"]:
  cert=c.get("exact_pf_certificate")
  if not cert:continue
  x={ast.literal_eval(s):int(v) for s,v in cert["positive_integer_vector"].items()}
  lo,hi,den=[int(cert[k]) for k in ("lower_numerator","upper_numerator","denominator")]
  assert min(x.values())>0
  certified_components.add(frozenset(x))
  for s in x:
   y=sum(v*x[t] for t,v in M[s].items() if t in x)
   assert lo*x[s]<=den*y<=hi*x[s]
  checked+=1
 assert certified_components==cyclic_components
 out["threshold"][d]={"states":len(M),"100_exact_count_terms_verified":True,"all_states_reachable":True,"all_cyclic_SCCs_certified":True,"SCC_PF_certificates_verified":checked,"first2":vals[:2]}
out["input_sha256"]={name:hashlib.sha256((base/name).read_bytes()).hexdigest() for name in ("exact_peptide_automata.json","exact_bridge_threshold_automata.json")}
(base/"automata_certificate_review.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps(out,indent=2))

