"""Independent exhaustive all-genotype check of shifted-register fibre theorem."""
import collections
import itertools
import json
from pathlib import Path
from independent_oracle import CODE, PAIRS, RY, CODONS, bfs, graph, register

def fib(k):
    a,b=0,1
    for _ in range(k): a,b=b,a+b
    return a

def fibre_formula(peptide, pairs, sample):
    n=len(sample)
    domains=[]
    for c,a in enumerate(peptide):
        cls=tuple(RY[b] for b in sample[3*c:3*c+3])
        codons=[codon for codon in CODONS[a] if tuple(RY[b] for b in codon)==cls]
        ds=[set(codon[j] for codon in codons) for j in range(3)]
        assert set(itertools.product(*ds))==set(map(tuple,codons)), (a,cls,'not product')
        domains.extend(ds)
    def dangerous(p): return 'A' if RY[sample[p]] else 'C'
    fixed_danger={p for p in range(n) if domains[p]=={dangerous(p)}}
    assert all(not(i in fixed_danger and j in fixed_danger) for i,j in pairs)
    free={p for p in range(n) if len(domains[p])==2}
    forced={j for i,j in pairs if i in fixed_danger}|{i for i,j in pairs if j in fixed_danger}
    free-=forced
    edges=set(tuple(sorted((i,j))) for i,j in pairs if i in free and j in free)
    adj={p:set() for p in free}
    for i,j in edges: adj[i].add(j);adj[j].add(i)
    unseen=set(free); lengths=[]
    while unseen:
        p=next(iter(unseen)); reached={p}; stack=[p]
        while stack:
            for q in adj[stack.pop()]-reached:
                reached.add(q); stack.append(q)
        unseen-=reached
        assert sum(len(adj[q]) for q in reached)//2==len(reached)-1
        assert max(map(lambda q:len(adj[q]),reached))<=2
        lengths.append(len(reached))
    return sorted(lengths), __import__('math').prod(fib(k+2) for k in lengths)

def main():
    n=6; summary={}; examples={}
    for d in range(1,n+1):
        p0=register(n,0); pd=register(n,d); pairs=tuple(sorted(set(p0+pd)))
        groups=collections.defaultdict(list); counts=collections.Counter()
        for types in itertools.product(sorted(PAIRS),repeat=n):
            x=['?']*(2*n)
            for (i,j),pair in zip(p0,types): x[i],x[j]=pair
            x=''.join(x)
            if any(x[i]+x[j] not in PAIRS for i,j in pd): continue
            peptide=tuple(CODE[x[i:i+3]] for i in range(0,len(x),3))
            if '*' in peptide: continue
            groups[(peptide,tuple(RY[c] for c in x))].append(x)
        for (peptide,classes),vs in groups.items():
            counts['nonempty_orientation_fibres']+=1
            lengths,size=fibre_formula(peptide,pairs,vs[0])
            assert size==len(vs),(d,peptide,classes,lengths,size,len(vs))
            if d==1: assert max(lengths,default=0)<=1
            if d==2: assert max(lengths,default=0)<=3
            counts['vertices']+=len(vs);counts['max_fibre_size']=max(counts['max_fibre_size'],len(vs))
            counts['max_free_path']=max(counts['max_free_path'],max(lengths,default=0))
            adj=graph(vs); ds,_=bfs(adj,0)
            assert len(ds)==len(vs)
            counts['edges']+=sum(map(len,adj))//2
            key=str(max(lengths,default=0))
            if key not in examples.setdefault(str(d),{}):
                examples[str(d)][key]={'peptide':''.join(peptide),'sequence':vs[0],'free_path_lengths':lengths,'fibre_size':size,'sequences':vs}
        counts['peptides_with_intersection']=len({p for p,cl in groups})
        summary[str(d)]=dict(counts)
    root=Path(__file__).resolve().parent
    (root/'register_fibonacci_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    (root/'register_fibonacci_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
