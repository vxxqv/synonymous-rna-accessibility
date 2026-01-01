"""Exact quotient and Hamming paths for overlapping STANDARD-code obligations.

Occurrences are (ordered three distinct positions, whole-codon complement).
Arbitrary triplets are mathematical constraints, not necessarily physical ORFs.
Pair obligations use underlying RNA bases, not complemented coding coordinates.
"""
from collections import deque
from pathlib import Path
import itertools
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'exact'))
from independent_oracle import CODE, CODONS, PAIRS

BASES = 'ACGU'
COMP = dict(zip('ACGU', 'UGCA'))


def read(sequence, occurrence):
    positions, complement = occurrence
    return ''.join(COMP[sequence[p]] if complement else sequence[p] for p in positions)


def validate(n, occurrences, peptide, edges, fixed):
    if not isinstance(n, int) or n < 0 or len(occurrences) != len(peptide):
        raise ValueError('nonnegative length and one sense amino acid per occurrence required')
    for (positions, complement), aa in zip(occurrences, peptide):
        if aa not in CODONS or len(positions) != 3 or len(set(positions)) != 3 or not isinstance(complement, bool):
            raise ValueError('distinct codon positions, Boolean complement and sense amino acid required')
        if any(not isinstance(p, int) or not 0 <= p < n for p in positions):
            raise ValueError('invalid codon coordinate')
    for i, j in edges:
        if i == j or not 0 <= i < n or not 0 <= j < n:
            raise ValueError('invalid pair edge')
    if any(not 0 <= p < n or b not in BASES for p, b in fixed.items()):
        raise ValueError('invalid fixed coordinate')


def valid(sequence, occurrences, peptide, edges=(), fixed=None):
    return (all(CODE.get(read(sequence, occ)) == aa for occ, aa in zip(occurrences, peptide))
            and all(sequence[i] + sequence[j] in PAIRS for i, j in edges)
            and all(sequence[p] == b for p, b in (fixed or {}).items()))


def sector(sequence, occurrences, peptide, edges):
    touched = {p for edge in edges for p in edge}
    paired = [p for p in range(len(sequence)) if p in touched]
    return (tuple(read(sequence, occ).startswith('AG') for occ, aa in zip(occurrences, peptide) if aa == 'S'),
            tuple(sequence[p] in 'AG' for p in paired))


def compile_sector(n, occurrences, peptide, edges, fixed, signature):
    """Exact coordinate quotient in a fixed Ser-family/paired-RY sector.

    A coordinate of size2 uses one Boolean variable for its individual bases.
    Size3/4 coordinates use one for R/Y; size1 coordinates use no variable.
    Returned bucket sets describe every original sequence in each fibre.
    """
    edges, fixed = tuple(edges), dict(fixed)
    validate(n, occurrences, peptide, edges, fixed)
    families, orientation = signature
    touched = {p for edge in edges for p in edge}
    paired = [p for p in range(n) if p in touched]
    if len(families) != sum(aa == 'S' for aa in peptide) or len(orientation) != len(paired):
        raise ValueError('complete Ser and paired-site signature required')
    domains = [set(BASES) for _ in range(n)]
    gates = []
    family_index = 0
    for occ, aa in zip(occurrences, peptide):
        positions, complement = occ
        codons = CODONS[aa]
        if aa == 'S':
            family = families[family_index]
            family_index += 1
            codons = tuple(c for c in codons if c.startswith('AG') == family)
        for local, p in enumerate(positions):
            allowed = {COMP[c[local]] if complement else c[local] for c in codons}
            domains[p].intersection_update(allowed)
        if aa in ('L', 'R'):
            trigger = 'U' if aa == 'L' else 'A'
            trigger = COMP[trigger] if complement else trigger
            gates.append((positions[0], trigger, positions[2], not complement))
    for p, purine in zip(paired, orientation):
        domains[p].intersection_update('AG' if purine else 'CU')
    for p, base in fixed.items():
        domains[p].intersection_update(base)
    result = {'domains': tuple(''.join(sorted(d)) for d in domains), 'buckets': [],
              'variable_sites': [], 'site_variables': {}, 'clauses': [], 'unsat': any(not d for d in domains)}
    if result['unsat']:
        return result
    for p, domain in enumerate(domains):
        bases = ''.join(sorted(domain))
        if len(domain) == 1:
            buckets = (bases,)
        elif len(domain) == 2:
            buckets = tuple(bases)
        else:
            buckets = (''.join(b for b in bases if b in 'CU'), ''.join(b for b in bases if b in 'AG'))
        result['buckets'].append(buckets)
        if len(buckets) == 2:
            result['site_variables'][p] = len(result['variable_sites'])
            result['variable_sites'].append(p)

    def predicate(p, test):
        truths = []
        for bucket in result['buckets'][p]:
            values = {bool(test(base)) for base in bucket}
            if len(values) != 1:
                raise AssertionError('predicate is not constant on a coordinate fibre')
            truths.append(values.pop())
        if len(set(truths)) == 1:
            return truths[0]
        variable = result['site_variables'][p] + 1
        return variable if truths[1] else -variable

    def negate(literal):
        return not literal if isinstance(literal, bool) else -literal

    def clause(*literals):
        if any(lit is True for lit in literals):
            return
        literals = set(lit for lit in literals if lit is not False)
        if any(-lit in literals for lit in literals):
            return
        if not literals:
            result['unsat'] = True
        else:
            result['clauses'].append(tuple(sorted(literals)))

    for p, trigger, q, desired_purine in gates:
        clause(negate(predicate(p, lambda b, trigger=trigger: b == trigger)),
               predicate(q, lambda b, desired_purine=desired_purine: (b in 'AG') == desired_purine))
    orientation_by_site = dict(zip(paired, orientation))
    for i, j in edges:
        if orientation_by_site[i] == orientation_by_site[j]:
            result['unsat'] = True
        clause(negate(predicate(i, lambda b: b in 'AC')),
               negate(predicate(j, lambda b: b in 'AC')))
    components, frozen = implication_components(result)
    result['unsat'] = result['unsat'] or any(
        components[2*v] == components[2*v+1]
        for v in range(len(result['variable_sites'])))
    result['frozen_variables'] = frozen
    return result


def project(compiled, sequence):
    return tuple(next(b for b, bucket in enumerate(compiled['buckets'][p]) if sequence[p] in bucket)
                 for p in compiled['variable_sites'])


def accepts(compiled, bits):
    return not compiled['unsat'] and all(any(bool(bits[abs(lit)-1]) == (lit > 0) for lit in clause)
                                       for clause in compiled['clauses'])


def implication_components(compiled):
    """Literal SCC labels and variables in SCCs with >=2 distinct variables."""
    n = len(compiled['variable_sites'])
    graph, reverse = [[] for _ in range(2*n)], [[] for _ in range(2*n)]
    def node(lit):
        return 2*(abs(lit)-1) + int(lit < 0)
    for clause in compiled['clauses']:
        a, b = node(clause[0]), node(clause[-1])
        for u, v in ((a ^ 1, b), (b ^ 1, a)):
            graph[u].append(v)
            reverse[v].append(u)
    seen, finish = set(), []
    for start in range(2*n):
        if start in seen:
            continue
        seen.add(start)
        stack = [(start, 0)]
        while stack:
            u, k = stack[-1]
            if k == len(graph[u]):
                stack.pop()
                finish.append(u)
            else:
                v = graph[u][k]
                stack[-1] = u, k+1
                if v not in seen:
                    seen.add(v)
                    stack.append((v, 0))
    component = [-1]*(2*n)
    members = []
    for start in reversed(finish):
        if component[start] != -1:
            continue
        label = len(members)
        component[start] = label
        stack, group = [start], []
        while stack:
            u = stack.pop()
            group.append(u)
            for v in reverse[u]:
                if component[v] == -1:
                    component[v] = label
                    stack.append(v)
        members.append(group)
    nontrivial = [len({u//2 for u in group}) > 1 for group in members]
    frozen = tuple(v for v in range(n) if nontrivial[component[2*v]] or nontrivial[component[2*v+1]])
    return tuple(component), frozen


def frozen_variables(compiled):
    """Variables in an implication SCC containing >=2 distinct variables."""
    if 'frozen_variables' in compiled:
        return compiled['frozen_variables']
    return implication_components(compiled)[1]


def endpoint_mutation_order(n, occurrences, peptide, edges, fixed, x, y):
    """Return coordinate flips in O(n+c+|E|) operations, or None if separated.

    The result records a compact nucleotide Hamming geodesic. It does not
    materialize all intermediate n-letter sequence strings.
    """
    edges, fixed = tuple(edges), dict(fixed)
    validate(n, occurrences, peptide, edges, fixed)
    if (len(x) != n or len(y) != n
            or any(base not in BASES for sequence in (x, y) for base in sequence)
            or not valid(x, occurrences, peptide, edges, fixed)
            or not valid(y, occurrences, peptide, edges, fixed)):
        raise ValueError('valid endpoints required')
    signature = sector(x, occurrences, peptide, edges)
    if signature != sector(y, occurrences, peptide, edges):
        return None
    compiled = compile_sector(n, occurrences, peptide, edges, fixed, signature)
    a, b = project(compiled, x), project(compiled, y)
    if any(a[v] != b[v] for v in frozen_variables(compiled)):
        return None
    changed = [v for v in range(len(a)) if a[v] != b[v]]
    changed_set = set(changed)
    relations = set()
    for clause in compiled['clauses']:
        local = sorted({abs(lit)-1 for lit in clause} & changed_set)
        for mask in range(1 << len(local)):
            values = {v: b[v] if mask & (1 << k) else a[v]
                      for k, v in enumerate(local)}
            if not any(bool(values.get(abs(lit)-1, a[abs(lit)-1])) == (lit > 0)
                       for lit in clause):
                if len(local) != 2 or mask not in (1, 2):
                    raise AssertionError('valid endpoints must give only binary event precedence')
                relations.add((local[0], local[1]) if mask == 2 else (local[1], local[0]))
    out = {v: [] for v in changed}
    indegree = {v: 0 for v in changed}
    for u, v in relations:
        out[u].append(v)
        indegree[v] += 1
    queue = deque(v for v in changed if not indegree[v])
    order = []
    while queue:
        v = queue.popleft()
        order.append(v)
        for w in out[v]:
            indegree[w] -= 1
            if not indegree[w]:
                queue.append(w)
    if len(order) != len(changed):
        raise AssertionError('SCC signature sufficiency failed')
    sites = [compiled['variable_sites'][v] for v in order]
    changed_sites = set(sites)
    for p in range(n):
        if x[p] != y[p] and p not in changed_sites:
            sites.append(p)
    return tuple(sites)


def endpoint_path(n, occurrences, peptide, edges, fixed, x, y):
    """Materialize a shortest path; writing all sequence states can cost O(n^2)."""
    order = endpoint_mutation_order(n, occurrences, peptide, edges, fixed, x, y)
    if order is None:
        return None
    sequence, path = list(x), [x]
    for p in order:
        sequence[p] = y[p]
        path.append(''.join(sequence))
    return path
