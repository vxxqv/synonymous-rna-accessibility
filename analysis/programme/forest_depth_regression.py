from pathlib import Path
from hashlib import sha256
import json
import sys
import time

root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'analysis/expansion'))
from graph_compiler import forest_extensions
n=1200
edges={(i,i+1) if i%2==0 else(i+1,i) for i in range(n-1)}
start=time.perf_counter()
result=forest_extensions(n,edges)
row=[1]
for m in range(1,n+1):
    old=row;row=[0]
    for k in range(1,m+1):row.append(row[-1]+old[m-k])
assert result==row[-1]
receipt={'status':'passed','events':n,'poset_height':2,'relations':n-1,
 'source_sha256':sha256((root/'analysis/expansion/graph_compiler.py').read_bytes()).hexdigest(),
 'independent_oracle':'Entringer triangle recurrence for alternating permutations',
 'exact_equal':True,'answer_bits':result.bit_length(),'recursion_limit':sys.getrecursionlimit(),
 'seconds':time.perf_counter()-start}
(root/'results/programme/forest_depth_regression.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
