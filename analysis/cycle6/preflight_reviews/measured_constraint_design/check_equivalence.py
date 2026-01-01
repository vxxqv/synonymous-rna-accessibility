"""Exact symbolic-row test of the supplied FAS factor type, no outcomes."""
if not __debug__:raise RuntimeError('Optimization refused')
from pathlib import Path
from itertools import product
from fractions import Fraction as Q
import hashlib,json
def require(x,msg):
 if not x:raise RuntimeError(msg)
def rank(a):
 a=[list(map(Q,r)) for r in a];row=0
 for c in range(len(a[0])):
  pivot=next((i for i in range(row,len(a)) if a[i][c]),None)
  if pivot is None:continue
  a[row],a[pivot]=a[pivot],a[row];p=a[row][c];a[row]=[x/p for x in a[row]]
  for i in range(len(a)):
   if i!=row and a[i][c]:
    p=a[i][c];a[i]=[x-p*y for x,y in zip(a[i],a[row])]
  row+=1
 return row
Xc=[];Xn=[]
for p,q,u,v,w in product(range(3),range(3),range(2),range(2),range(2)):
 a,b=[(0,0),(0,1),(1,1)][p];c,d=[(1,1),(0,1),(0,0)][q]
 require((int(p==1),int(p==2))==(b-a,a),'First Leu forward identity')
 require((int(q==1),int(q==2))==(d-c,1-d),'Second Leu forward identity')
 require((a,b)==(int(p==2),int(p==1)+int(p==2)),'First inverse')
 require((c,d)==(1-int(q==1)-int(q==2),1-int(q==2)),'Second inverse')
 require(1+a-b==[1,0,1][p],'Nucleotide-additive gate valley')
 Xc.append([1,int(p==1),int(p==2),int(q==1),int(q==2),u,v,w]);Xn.append([1,a,b,c,d,u,v,w])
ranks=[rank(Xc),rank(Xn),rank([a+b for a,b in zip(Xc,Xn)])]
require(ranks==[8,8,8],'Additive span disagreement')
out={'status':'EXACT_DOMAIN_SPECIFIC_ADDITIVE_SPANS_EQUAL','factor_rows':len(Xc),'ranks_codon_nucleotide_joined':ranks,'explicit_forward_and_inverse_identities_checked':True,'P3_gate_valley_is_also_nucleotide_additive':True,'biological_outcomes_used':False,'independent_source_sequence_reconstruction':False,'scope':'Algebra on source-lane supplied factor states; no new scientific result or generic all-domains equivalence','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
p=Path(__file__).with_name('baseline_equivalence_checks.json')
if p.exists():raise RuntimeError('No overwrite')
p.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps(out,indent=2))
