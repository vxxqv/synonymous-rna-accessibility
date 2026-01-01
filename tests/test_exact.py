import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from accessibility import *
import pytest

def bfs(graph, root):
    distance, paths, queue = {root:0}, {root:1}, [root]
    for u in queue:
        for v in graph[u]:
            if v not in distance:
                distance[v]=distance[u]+1;paths[v]=0;queue.append(v)
            if distance[v]==distance[u]+1:
                paths[v]+=paths[u]
    return distance,paths

def test_unique_three_mutation_route():
    seq=enumerate_space('LG',[(0,5)])
    graph=adjacency(seq)
    assert len(seq)==8 and sum(map(len,graph))//2==12
    result=endpoint_certificate('CUCGGG','UUAGGA',[(0,5)])
    assert result['shortest_paths']=='1'
    assert result['chains']==[[2,0,5]]

def test_single_codon_and_pairing_oracles():
    for aa in CODONS:
        for pairs in [(),((0,1),),((0,2),),((1,2),)]:
            seq=enumerate_space(aa,pairs)
            graph=adjacency(seq)
            for i,x in enumerate(seq):
                distances,paths=bfs(graph,i)
                for j,y in enumerate(seq):
                    cert=endpoint_certificate(x,y,pairs)
                    assert cert['connected']==(j in distances)
                    if j in distances:
                        assert cert['distance']==distances[j]
                        assert int(cert['shortest_paths'])==paths[j]

def test_same_endpoint_fragmentation():
    seq=enumerate_space('LG',[(0,5)])
    graph=adjacency(seq)
    all_keep=range(len(seq))
    assert excess_fragmentation(seq,[(0,5)],graph,all_keep)['excess']==0
    assert excess_fragmentation(seq,[(0,5)],graph,[])['excess'] is None
    # With two endpoints disconnected by deleting intermediates, baseline
    # connectivity remains one while filtered connectivity is one half.
    selected=[seq.index(x) for x in ('CUCGGG','UUAGGA')]
    assert excess_fragmentation(seq,[(0,5)],graph,selected)['excess']==0.5

def test_fixed_coordinates():
    seq=enumerate_space('LG',[(0,5)],fixed={2:'A'})
    assert seq and all(s[2]=='A' for s in seq)
    for x in seq:
        for y in seq:
            result=endpoint_certificate(x,y,[(0,5)])
            if result['connected']:
                assert all(z[2]=='A' for z in result['witness'])

def test_invalid_inputs():
    with pytest.raises(ValueError):matching([(0,2),(1,2)],3)
    with pytest.raises(MemoryError):enumerate_space('LLLLLLLL',cap=10)
    with pytest.raises(ValueError):endpoint_certificate('AUG','GGG',[])

def test_register_fixture():
    for d,expected in [(1,{0:24,1:16,2:6,4:6}),(2,{0:52,1:32})]:
        hist=Counter()
        for a,b in product(CODONS,repeat=2):
            c0=enumerate_space(a+b,register_pairs(3,0))
            cd=enumerate_space(a+b,register_pairs(3,d))
            delta=cross_layer_distance(c0,cd)
            if delta is not None:hist[delta]+=1
        assert dict(hist)==expected
