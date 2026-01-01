"""Independent invented fixtures and exhaustive simple-path reference."""
import sys
if not __debug__:
    raise RuntimeError('Optimized Python is not allowed')
import copy
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import itertools
import json
import math
import subprocess
import tempfile
import verify_result as v


def rejected(function):
    try:
        function()
    except (AssertionError, ValueError):
        return
    raise AssertionError('Expected rejection')


def paths_oracle(scores, adj, a, b):
    def search(path):
        if path[-1] == b:
            return min(scores[k] for k in path)
        return max([search(path+[j]) for j in adj[path[-1]] if j not in path], default=-math.inf)
    return search([a])


def main():
    feasible={'status':'feasible','forced_orders':1,'exact_checked':True,
              'beta':['0','1'],'latent':['0','1'],'response':['0','1']}
    assert v.certificate(feasible,[[1,0],[1,1]],[Q(0),Q(1)],[Q(0),Q(1)])['status']=='feasible'
    constant={'status':'feasible','forced_orders':0,'exact_checked':True,
              'beta':['0'],'latent':['0','0'],'response':['2','2']}
    v.certificate(constant,[[1],[1]],[Q(0),Q(2)],[None,Q(3)])
    # A reciprocal-sign Q2 has two forced orders with opposite additive differences.
    X=[[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    dual={'status':'infeasible','forced_orders':4,'exact_checked':True,'dual':[
        {'lower_vertex':1,'higher_vertex':0,'weight':'1/2','difference':[0,0,-1]},
        {'lower_vertex':2,'higher_vertex':3,'weight':'1/2','difference':[0,0,1]}]}
    L=[Q(1),Q(0),Q(0),Q(1)]
    assert v.certificate(dual,X,L,L)['support_size']==2
    for key in ('weight','difference','higher_vertex'):
        bad=copy.deepcopy(dual)
        bad['dual'][0][key]={'weight':'1/3','difference':[0,1,-1],'higher_vertex':3}[key]
        rejected(lambda:v.certificate(bad,X,L,L))
    # Touching intervals are not forced, and ties constrain equal response.
    rejected(lambda:v.certificate(feasible,[[1,0],[1,1]],[Q(0),Q(1)],[Q(1),Q(2)]))
    bad=copy.deepcopy(constant);bad['response']=['1','2']
    rejected(lambda:v.certificate(bad,[[1],[1]],[Q(0),Q(0)],[Q(3),Q(3)]))
    graph_cases=pair_values=0
    # Every undirected graph on <=4 vertices, four value assignments per vertex,
    # including rational negatives, equal levels, and unknown optimistic infinity.
    for n in range(1,5):
        possible=list(itertools.combinations(range(n),2))
        for mask in range(1<<len(possible)):
            edges=[e for k,e in enumerate(possible) if mask>>k&1]
            adj=v.adjacency(n,edges)
            for scores in itertools.product([Q(-1,2),Q(0),Q(1),math.inf],repeat=n):
                actual=v.capacities(scores,adj)
                graph_cases+=1
                for a in range(n):
                    for b in range(a,n):
                        assert actual[a][b]==paths_oracle(scores,adj,a,b)
                        if a!=b:
                            v.witness(a,b,scores,adj,actual[a][b])
                        pair_values+=1
    # Unknown center permits an optimistic route, weak equality retains a path.
    adj=v.adjacency(3,[(0,1),(1,2)])
    assert v.capacities([Q(1),math.inf,Q(1)],adj)[0][2]==1
    assert v.capacities([Q(1),Q(0),Q(1)],adj)[0][2]==0
    assert 2 in v.component(adj,[Q(1),Q(0),Q(1)],0,0)
    assert 2 not in v.component(adj,[Q(1),Q(0),Q(1)],0,0,strict=True)
    # Optimization must fail at module start, even absent arguments.
    guards=[]
    for flag in ('-O','-OO'):
        p=subprocess.run([sys.executable,flag,str(Path(v.__file__).resolve())],capture_output=True,text=True)
        assert p.returncode!=0 and 'Optimized Python is not allowed' in p.stderr
        guards.append(flag)
    receipt={'status':'PASS','graph_score_cases':graph_cases,'simple_path_capacity_values':pair_values,
             'certificate_fixtures':'feasible, plateau/infinite upper, exact dual, altered weights/differences/orders, weak touching intervals',
             'optimization_rejected':guards,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in [Path(__file__),Path(v.__file__)]},
             'scope':'Synthetic finite oracle checks, not a theorem or source-count authentication.'}
    out=Path(__file__).parent/'synthetic_receipt.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
