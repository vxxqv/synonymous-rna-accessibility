"""Exact standard-code compatibility and endpoint certificates.

All nucleotide positions are zero based. The model permits one matching,
one standard reading frame, and fixed coordinates. It does not model fitness.
"""
from collections import Counter, deque
from itertools import product
from math import factorial

BASES = 'UCAG'
AA = ('FFLLSSSSYY**CC*W'
      'LLLLPPPPHHQQRRRR'
      'IIIMTTTTNNKKSSRR'
      'VVVVAAAADDEEGGGG')
CODE = dict(zip((''.join(c) for c in product(BASES, repeat=3)), AA))
CODONS = {a: tuple(c for c in CODE if CODE[c] == a) for a in sorted(set(AA)-{'*'})}
PAIRS = frozenset(('AU','UA','GC','CG','GU','UG'))

def rna(sequence):
    sequence = sequence.upper().replace('T','U')
    if set(sequence) - set(BASES):
        raise ValueError('Only unambiguous A C G U/T bases are accepted')
    return sequence

def translate(sequence):
    sequence = rna(sequence)
    if len(sequence) % 3:
        raise ValueError('Incomplete codon')
    return ''.join(CODE[sequence[i:i+3]] for i in range(0,len(sequence),3))

def matching(pairs, length):
    pairs = tuple(sorted((min(i,j), max(i,j)) for i,j in pairs))
    flat = [v for edge in pairs for v in edge]
    if len(set(flat)) != len(flat) or any(i < 0 or i >= length for i in flat) or any(i == j for i,j in pairs):
        raise ValueError('Pair table must be a matching on the sequence')
    return pairs

def compatible(sequence, pairs):
    return all(sequence[i]+sequence[j] in PAIRS for i,j in pairs)

def signature(sequence, pairs):
    return tuple(sequence[i] in 'AG' for i,j in pairs) + tuple(
        sequence[i:i+2] == 'AG' for i in range(0,len(sequence),3)
        if CODE[sequence[i:i+3]] == 'S')

def enumerate_space(peptide, pairs=(), fixed=None, cap=50000):
    if any(a not in CODONS for a in peptide):
        raise ValueError('Peptide must use standard sense amino acids')
    pairs = matching(pairs,3*len(peptide))
    fixed = fixed or {}
    if any(i < 0 or i >= 3*len(peptide) or b not in BASES for i,b in fixed.items()):
        raise ValueError('Invalid fixed coordinate')
    domains = [tuple(c for c in CODONS[a] if all(3*k+j not in fixed or fixed[3*k+j] == c[j] for j in range(3))) for k,a in enumerate(peptide)]
    size = 1
    for domain in domains:
        size *= len(domain)
    if size > cap:
        raise MemoryError(f'Cartesian domain {size} exceeds cap {cap}')
    return [''.join(c) for c in product(*domains) if compatible(''.join(c),pairs)]

def adjacency(sequences):
    sequences = tuple(sequences)
    if len(set(sequences)) != len(sequences):
        raise ValueError('Duplicate graph vertices')
    lookup = {s:i for i,s in enumerate(sequences)}
    graph = [[] for _ in sequences]
    for i,s in enumerate(sequences):
        for k,old in enumerate(s):
            for new in BASES:
                if new != old:
                    j = lookup.get(s[:k]+new+s[k+1:])
                    if j is not None:
                        graph[i].append(j)
    return graph

def components(graph, keep=None):
    keep = set(range(len(graph))) if keep is None else set(keep)
    labels, sizes = {}, []
    for root in sorted(keep):
        if root in labels:
            continue
        index, size, queue = len(sizes), 0, [root]
        labels[root] = index
        for node in queue:
            size += 1
            for other in graph[node]:
                if other in keep and other not in labels:
                    labels[other] = index
                    queue.append(other)
        sizes.append(size)
    return labels, sizes

def connected_probability(sizes):
    n = sum(sizes)
    return None if n == 0 else sum(s*s for s in sizes)/(n*n)

def excess_fragmentation(sequences, pairs, graph, keep):
    keep = tuple(keep)
    if not keep:
        return {'n':0,'a_filtered':None,'a_compatibility_same_endpoints':None,'excess':None}
    baseline = Counter(signature(sequences[i],pairs) for i in keep)
    _, sizes = components(graph,keep)
    numerator = sum(n*n for n in sizes)
    denominator = sum(n*n for n in baseline.values())
    return {'n':len(keep),'components':len(sizes),'sizes':sorted(sizes,reverse=True),
            'same_component_pairs_filtered':str(numerator),'same_component_pairs_compatibility':str(denominator),
            'a_filtered':numerator/len(keep)**2,
            'a_compatibility_same_endpoints':denominator/len(keep)**2,
            'excess':1-numerator/denominator}

def endpoint_certificate(start, end, pairs):
    start, end = rna(start), rna(end)
    if len(start) != len(end) or translate(start) != translate(end) or '*' in translate(start):
        raise ValueError('Endpoints must encode the same sense peptide')
    pairs = matching(pairs,len(start))
    if not compatible(start,pairs) or not compatible(end,pairs):
        raise ValueError('Endpoints violate pairing constraints')
    if signature(start,pairs) != signature(end,pairs):
        return {'connected':False,'shortest_paths':'0'}
    changed = [i for i in range(len(start)) if start[i] != end[i]]
    precedence = set()
    factors = [tuple(range(i,i+3)) for i in range(0,len(start),3)] + list(pairs)
    for factor in factors:
        local = [i for i in factor if i in changed]
        if len(local) < 2:
            continue
        valid_subsets = []
        for bits in product((0,1),repeat=len(local)):
            x = list(start)
            switched = {i for i,b in zip(local,bits) if b}
            for i in switched:
                x[i] = end[i]
            if len(factor) == 3:
                ok = CODE[''.join(x[i] for i in factor)] == CODE[start[factor[0]:factor[0]+3]]
            else:
                ok = x[factor[0]]+x[factor[1]] in PAIRS
            if ok:
                valid_subsets.append(switched)
        for i in local:
            for j in local:
                if i != j and all(j not in subset or i in subset for subset in valid_subsets):
                    precedence.add((i,j))
    outgoing, incoming = {i:[] for i in changed}, {i:[] for i in changed}
    for i,j in precedence:
        outgoing[i].append(j)
        incoming[j].append(i)
    if any(len(outgoing[i]) > 1 or len(incoming[i]) > 1 for i in changed):
        raise AssertionError('Non-chain precedence in the stated model')
    chains, visited = [], set()
    for root in changed:
        if incoming[root]:
            continue
        chain, node = [], root
        while True:
            if node in visited:
                raise AssertionError('Precedence cycle')
            chain.append(node)
            visited.add(node)
            if not outgoing[node]:
                break
            node = outgoing[node][0]
        chains.append(chain)
    if len(visited) != len(changed) or any(len(c)>3 for c in chains):
        raise AssertionError('Cycle or chain length exceeds theorem bound')
    count = factorial(len(changed))
    for chain in chains:
        count //= factorial(len(chain))
    # Every linear extension gives a geodesic. Concatenating chains supplies one.
    witness, x = [start], list(start)
    for chain in chains:
        for node in chain:
            x[node] = end[node]
            y = ''.join(x)
            if translate(y) != translate(start) or not compatible(y,pairs):
                raise AssertionError('Invalid witness')
            witness.append(y)
    return {'connected':True,'distance':len(changed),'chains':chains,
            'precedence':sorted(precedence),'shortest_paths':str(count),'witness':witness}

def register_pairs(n, shift=0):
    """Two arms concatenated without the unchanged loop, 0 based."""
    if n < 1 or not 0 <= shift <= n:
        raise ValueError('Require 0 <= shift <= arm length')
    return tuple((i, n+(n-1+shift-i)) for i in range(shift,n))

def cross_layer_distance(layer0, layer1):
    if not layer0 or not layer1:
        return None
    return min(sum(a!=b for a,b in zip(x,y)) for x in layer0 for y in layer1)
