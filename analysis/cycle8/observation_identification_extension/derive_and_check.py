"""Exact synthetic checks and complete finite observation-subset inventory."""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import ast,hashlib,json,re
if not __debug__:raise RuntimeError('Optimized execution refused')
P=Path(__file__).resolve().parent
N,C,PAMOUNT=F(2),F(3),F(5)
GROUPS={'A':('s','a'),'E':('e','d_C','J'),'D':('d_N',),'F':('f',),'B':('k','d_P','H')}
OBS=tuple(v for g in GROUPS.values() for v in g)
CLASS={v:g for g,vs in GROUPS.items() for v in vs}
def physical(x,y,z):return {'s':N*(x+y),'e':x,'d_N':y,'d_C':N*x/C,'k':PAMOUNT*z/C,'d_P':z}
def observations(x,y,z):return {**physical(x,y,z),'a':x+y,'J':N*x,'H':PAMOUNT*z,'f':x/(x+y)}
def base_value(name,v):
    return v/N if name in ('s','J') else v*C/N if name=='d_C' else v*C/PAMOUNT if name=='k' else v/PAMOUNT if name=='H' else v
def recover(sub,values):
    byclass={}
    for name in sub:
        g=CLASS[name];v=base_value(name,values[name])
        if (g=='F' and not F(0)<v<F(1)) or (g!='F' and v<=0):raise ValueError('Observation outside declared positive domain')
        if g in byclass and byclass[g]!=v:raise ValueError('Inconsistent redundant observation')
        byclass[g]=v
    if 'B' not in byclass or len(set(byclass)-{'B'})<2:raise ValueError('Insufficient')
    g1,g2=[g for g in ('A','E','D','F') if g in byclass][:2]
    u,v=byclass[g1],byclass[g2]
    if (g1,g2)==('A','E'):x,y=v,u-v
    elif (g1,g2)==('A','D'):x,y=u-v,v
    elif (g1,g2)==('A','F'):x,y=u*v,u*(1-v)
    elif (g1,g2)==('E','D'):x,y=u,v
    elif (g1,g2)==('E','F'):x,y=u,u*(1-v)/v
    elif (g1,g2)==('D','F'):x,y=u*v/(1-v),u
    else:raise AssertionError('Missing nuclear class pair')
    z=byclass['B']
    if min(x,y,z)<=0:raise ValueError('Outside strict positive domain')
    reconstructed=observations(x,y,z)
    if any(reconstructed[n]!=values[n] for n in sub):raise ValueError('Overidentified inconsistency')
    return x,y,z
def sufficient(sub):return bool(set(sub)&set(GROUPS['B'])) and len({CLASS[x] for x in sub}-{'B'})>=2
def alternate(sub,x,y,z):
    cs={CLASS[n] for n in sub}
    if 'B' not in cs:return x,y,z+1
    cs-= {'B'}
    if not cs:return x+1,y,z
    if cs=={'A'}:return x/2,y+x/2,z
    if cs=={'E'}:return x,y+1,z
    if cs=={'D'}:return x+1,y,z
    if cs=={'F'}:return 2*x,2*y,z
    raise AssertionError('Not an insufficient set')
subsets=[];minimal=[];checks=0
anchors=[(F(2),F(3),F(5)),(F(1,7),F(11,3),F(2,9)),(F(3),F(1,10**6),F(4)),(F(1,10**6),F(3),F(4)),(F(7,2),F(7,2),F(1,11))]
for mask in range(1<<len(OBS)):
    sub=tuple(OBS[i] for i in range(len(OBS)) if mask>>i&1)
    good=sufficient(sub)
    for x,y,z in anchors:
        values=observations(x,y,z)
        if good:
            if recover(sub,values)!=(x,y,z):raise AssertionError(('inverse',sub))
        else:
            xx,yy,zz=alternate(sub,x,y,z)
            if min(xx,yy,zz)<=0 or (xx,yy,zz)==(x,y,z) or any(observations(xx,yy,zz)[n]!=values[n] for n in sub):raise AssertionError(('counterexample',sub))
        checks+=1
    isminimal=good and not any(sufficient(tuple(x for x in sub if x!=drop)) for drop in sub)
    row={'observations':sub,'globally_sufficient':good,'minimal':isminimal,'classes':sorted({CLASS[n] for n in sub})}
    subsets.append(row)
    if isminimal:minimal.append(row)
assert len(subsets)==1024 and sum(x['globally_sufficient'] for x in subsets)==805 and len(minimal)==51
assert all(len(x['observations'])==3 for x in minimal)

# Boundary and inconsistent-input checks, not extensions of the positive theorem.
for sub,values in [(('a','e','d_P'),{'a':F(2),'e':F(2),'d_P':F(3)}),(('s','a','e','d_P'),{'s':F(9),'a':F(5),'e':F(2),'d_P':F(3)}),(('a','e','f','d_P'),{'a':F(5),'e':F(2),'f':F(1,2),'d_P':F(3)})]:
    try:recover(sub,values)
    except ValueError:checks+=1
    else:raise AssertionError('Inconsistent/boundary input accepted')
for invalid_f in (F(-1),F(0),F(1),F(2)):
    try:recover(('e','f','d_P'),{'e':F(2),'f':invalid_f,'d_P':F(3)})
    except ValueError:checks+=1
    else:raise AssertionError('Invalid export fraction accepted')

# Captured-category matrix: full-length F, two latent paths within the SAME
# reported last-junction bucket, and invisible/unretained mass.
A=[[F(1),F(0),F(0),F(0)],[F(0),F(1),F(1),F(0)]]
def mv(m,v):return [sum(a*b for a,b in zip(row,v)) for row in m]
null=[(F(0),F(1),F(-1),F(0)),(F(0),F(0),F(0),F(1))]
assert all(mv(A,v)==[0,0] for v in null)
u=list(map(F,(2,3,1,1)));v=list(map(F,(2,1,3,1)))
assert mv(A,u)==mv(A,v)==[2,4] and sum(u)==sum(v)==7
# Hypothetical coding contrast assigns the full-length class and first path
# to one product, second path/unretained class to another. Assignment is NOT
# authenticated for the actual source.
coding=(F(1),F(1),F(0),F(0))
assert sum(a*b for a,b in zip(coding,u))==5 and sum(a*b for a,b in zip(coding,v))==3

# Demonstrate the source parser's last-junction information loss without
# generating nucleotide sequences or asserting these paths occur in the data.
tree=ast.parse(Path('work/research_cycle8/dao_gate/source_junction_parser.py').read_text())
fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='categorize_read')
ns={'re':re};exec(compile(ast.Module(body=[fn],type_ignores=[]),'source_function_only','exec'),ns)
class PredicateReference:
    def __getitem__(self,key):
        pair=(key.start,key.stop)
        if pair in ((10,12),(100,102)):return 'GT'
        if pair in ((38,40),(128,130)):return 'AG'
        return 'lastD' if pair==(94,100) else 'lastA' if pair==(130,136) else 'other'
same1=ns['categorize_read']('100M30D50M',30,PredicateReference())
same2=ns['categorize_read']('10M30D60M30D50M',60,PredicateReference())
assert same1==same2=='AG_junction_low_mismatch-lastD_100_lastA_130'

# Six eligible barcodes, each100reads and positive full-length count. The
# coordinatewise median fractions sum to3/10 rather than1.
barcode_counts=[(80,10,10),(80,10,10),(10,80,10),(10,80,10),(10,10,80),(10,10,80)]
medians=[]
for col in zip(*barcode_counts):
    vals=sorted(F(x,100) for x in col);medians.append((vals[2]+vals[3])/2)
assert medians==[F(1,10)]*3 and sum(medians)==F(3,10)
# Missing category in one preparation is absent, not a zero row before pivot:
# mean-present S=.1 versus actual preparationmean S=.05.
missing_prep={'source_pivot_F':F(9,10)/2+F(1)/2,'source_pivot_S':F(1,10),'correct_zero_completed_prep_mean_S':F(1,20)}
assert missing_prep['source_pivot_F']+missing_prep['source_pivot_S']==F(21,20)

# One idealized occupancy cannot replace a protein clock if residence is free.
for z in (F(1,2),F(1),F(5)):
    occupancy=F(7);tau=occupancy/(PAMOUNT*z)
    assert tau>0 and PAMOUNT*z*tau==occupancy

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
outputs={'groups':GROUPS,'observation_order':OBS,'total_subsets':1024,'sufficient_subsets':805,'insufficient_subsets':219,'minimal_sets':minimal,'subsets':subsets}
(P/'observation_subsets.json').write_text(json.dumps(serial(outputs),indent=2)+'\n')
examples={'category_matrix':A,'rank':2,'nullspace_basis':null,'mixtures':[u,v],'same_observations':[2,4],'same_total':7,'hypothetical_coding_contrast':coding,'coding_contrast_values':[5,3],'coding_assignment_status':'admissible abstract unresolved assignment, not authenticated source transcript mixture','parser_categories':[same1,same2],'parser_fixture_status':'symbolic predicate reference, not a nucleotide sequence or assertion that actual GT-AG paths occur','six_barcode_medians':medians,'sum_median_fractions':sum(medians),'missing_preparation_example':missing_prep}
(P/'rational_examples.json').write_text(json.dumps(serial(examples),indent=2)+'\n')
receipt={'status':'PASS_EXACT_SYNTHETIC_CHECKS_PENDING_FRESH_REVIEW','finite_assertion':'All1024subsets classified, every sufficient set reconstructed and every insufficient set given positive indistinguishable alternatives at five rational anchors. Universal proof separate.','subset_anchor_checks':checks,'normalization':'n=2,c=3,p=5 for finite tests; proof covers every fixed positive triple','new_nucleotide_sequences_generated':0,'source_phenotype_rows_read':0,'author_code_execution':'Only original categorize_read function AST on symbolic CIGAR/predicate reference, no script top-level or file I/O','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({**receipt,'minimal_set_count':len(minimal)},indent=2))
