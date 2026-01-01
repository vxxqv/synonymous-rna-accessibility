"""Exact six-state standard-code counterexamples, with no biological measurements."""
from collections import deque
from pathlib import Path
import hashlib, json

VERTICES = ('CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG')
EDGES = tuple((u,v) for i,u in enumerate(VERTICES) for v in VERTICES[i+1:]
              if sum(a != b for a,b in zip(u,v)) == 1)

def route(scores, start, end, mode, threshold=0):
    adjacency={v:[] for v in VERTICES}
    for a,b in EDGES:
        adjacency[a].append(b); adjacency[b].append(a)
    parent={start:None}; queue=deque([start])
    if scores[start] < threshold: return None
    while queue:
        u=queue.popleft()
        if u==end:
            result=[]
            while u is not None: result.append(u); u=parent[u]
            return result[::-1]
        for v in adjacency[u]:
            allowed=(mode=='threshold' or
                     (mode=='weak' and scores[v]>=scores[u]) or
                     (mode=='strict' and scores[v]>scores[u]))
            if allowed and scores[v]>=threshold and v not in parent:
                parent[v]=u;queue.append(v)
    return None

def main():
    if not __debug__: raise RuntimeError('verification requires assertions')
    from Bio.Data.CodonTable import standard_rna_table
    assert tuple(sorted(c for c,a in standard_rna_table.forward_table.items() if a=='L'))==VERTICES
    rows=[]
    for kind,middle in (('threshold_without_weak_ascent',2),('weak_without_strict_ascent',0)):
        scores={v:2 for v in VERTICES}
        scores.update(UUA=0,CUC=1,CUA=middle)
        paths={mode:route(scores,'UUA','CUC',mode) for mode in ('threshold','weak','strict')}
        assert paths['threshold'] is not None
        assert (paths['weak'] is not None)==(middle==0)
        assert paths['strict'] is None
        rows.append(dict(case=kind,scores=scores,paths=paths))
    result=dict(scope='Abstract complete six-state leucine graph; arbitrary illustrative scores, not measurements',
                vertices=VERTICES,edges=EDGES,cases=rows,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                prior='Schmiegelt and Krug 2023 DOI 10.1007/s00285-023-01882-z treats monotone accessibility; no novelty claim for graph orientation')
    output=Path(__file__).resolve().parents[2]/'results/programme/accessibility_concepts.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('Two complete-code counterexamples verified against the standard genetic code.')

if __name__=='__main__': main()
