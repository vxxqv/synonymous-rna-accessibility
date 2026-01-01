"""Complete codon-CSP enumeration with separate node/output resource gates.

Unlike Cartesian enumeration, pruning can recover small feasible sets inside
enormous synonymous products. A resource failure never returns a partial set.
"""
from accessibility import CODONS,PAIRS

def enumerate_context(peptide,context,variable_start,pairs,output_cap=4096,node_cap=200000):
    length=3*len(peptide);v=variable_start
    if v<0 or v+length>len(context):raise ValueError('Variable interval outside context')
    if any(a not in CODONS for a in peptide):raise ValueError('Standard sense peptide required')
    if any(a<0 or b>=len(context) or a>=b for a,b in pairs):raise ValueError('Invalid pair coordinates')
    unary={k:[] for k in range(len(peptide))};binary={};constant=True
    for a,b in pairs:
        left=v<=a<v+length;right=v<=b<v+length
        if not left and not right:constant=constant and context[a]+context[b] in PAIRS
        elif left and right:
            ca,cb=(a-v)//3,(b-v)//3;ia,ib=(a-v)%3,(b-v)%3
            if ca==cb:unary[ca].append((ia,ib,None))
            else:binary.setdefault((ca,cb),[]).append((ia,ib))
        else:
            x,y=(a,b) if left else (b,a);unary[(x-v)//3].append(((x-v)%3,None,context[y]))
    if not constant:return [],{'nodes':0,'largest_output':0,'status':'exact'}
    domains=[]
    for k,aa in enumerate(peptide):
        d=[c for c in CODONS[aa] if all((c[i]+c[j] if j is not None else c[i]+base) in PAIRS for i,j,base in unary[k])];domains.append(d)
    if any(not d for d in domains):return [],{'nodes':0,'largest_output':0,'status':'exact'}
    neighbors={k:set() for k in range(len(peptide))}
    for a,b in binary:neighbors[a].add(b);neighbors[b].add(a)
    order=sorted(range(len(peptide)),key=lambda k:(len(domains[k]),-len(neighbors[k]),k))
    assigned={};output=[];nodes=0
    def visit(depth):
        nonlocal nodes
        if depth==len(order):
            if len(output)>=output_cap:raise MemoryError(f'Compatible output exceeds {output_cap}; no partial answer returned')
            output.append(''.join(assigned[k] for k in range(len(peptide))));return
        k=order[depth]
        for c in domains[k]:
            nodes+=1
            if nodes>node_cap:raise MemoryError(f'Search nodes exceed {node_cap}; no partial answer returned')
            assigned[k]=c
            if all(all(assigned[a][i]+assigned[b][j] in PAIRS for i,j in edges) for (a,b),edges in binary.items() if a in assigned and b in assigned):visit(depth+1)
            del assigned[k]
    visit(0)
    return sorted(output),{'nodes':nodes,'largest_output':len(output),'status':'exact','order':order,'domain_sizes':[len(d) for d in domains]}
