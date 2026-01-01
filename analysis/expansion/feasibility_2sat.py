"""Linear-size exact standard-code feasibility compiler for arbitrary pair graphs.

The two bits per nucleotide are (purine, dangerous): A=11,G=10,C=01,U=00.
Boolean-bit flips are NOT the nucleotide mutation graph. This module addresses
existence and construction only. No auxiliary Boolean variables are required.
"""
import itertools
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'exact'))
from independent_oracle import CODONS

BITS = {'A': (1, 1), 'G': (1, 0), 'C': (0, 1), 'U': (0, 0)}
BASE = {value: key for key, value in BITS.items()}


def encode(sequence):
    return tuple(bit for base in sequence for bit in BITS[base])


def evaluate(clauses, bits):
    return all(any(bool(bits[abs(lit) - 1]) == (lit > 0) for lit in clause)
               for clause in clauses)


def _codon_clauses(codons):
    """All entailed unit clauses and clauses on two distinct variables."""
    allowed = tuple(encode(codon) for codon in codons)
    candidates = [(sign * i,) for i in range(1, 7) for sign in (-1, 1)]
    candidates += [(sign * i, other * j)
                   for i in range(1, 7) for j in range(i + 1, 7)
                   for sign in (-1, 1) for other in (-1, 1)]
    return tuple(clause for clause in candidates
                 if all(evaluate((clause,), bits) for bits in allowed))


CODON_CNF = {aa: _codon_clauses(codons) for aa, codons in CODONS.items()}


def compile_instance(peptide, edges, fixed=None):
    """Return (6*m, clauses) with signed, one-based DIMACS literals.

    Input positions are zero based; fixed positions map to RNA bases.
    Size <=66*m + 3*|H| + 2*|fixed| clauses, each of length at most two.
    """
    if any(aa not in CODONS for aa in peptide):
        raise ValueError('standard sense peptide required')
    edges = tuple(edges)
    fixed = {} if fixed is None else dict(fixed)
    n = 3 * len(peptide)
    for edge in edges:
        if len(edge) != 2:
            raise ValueError('pair must have two positions')
        i, j = edge
        if not isinstance(i, int) or not isinstance(j, int) or i == j or not (0 <= i < n and 0 <= j < n):
            raise ValueError('invalid pairing edge')
    if any(not isinstance(p, int) or not 0 <= p < n or b not in BITS
           for p, b in fixed.items()):
        raise ValueError('fixed positions must be valid RNA bases at valid sites')
    clauses = []
    for c, aa in enumerate(peptide):
        for clause in CODON_CNF[aa]:
            clauses.append(tuple((abs(lit) + 6 * c) * (1 if lit > 0 else -1)
                                 for lit in clause))
    for i, j in edges:
        ri, rj, di, dj = 2 * i + 1, 2 * j + 1, 2 * i + 2, 2 * j + 2
        clauses.extend(((ri, rj), (-ri, -rj), (-di, -dj)))
    for p, base in fixed.items():
        for b, bit in enumerate(BITS[base]):
            variable = 2 * p + b + 1
            clauses.append((variable if bit else -variable,))
    return 2 * n, tuple(clauses)


def solve_cnf(nvariables, clauses):
    """Iterative Kosaraju SCC solver; None means unsatisfiable, () is SAT."""
    graph = [[] for _ in range(2 * nvariables)]
    reverse = [[] for _ in graph]

    def node(literal):
        return 2 * (abs(literal) - 1) + int(literal < 0)

    for clause in clauses:
        if not 1 <= len(clause) <= 2:
            raise ValueError('one or two literals per clause required')
        if any(not isinstance(lit, int) or not 1 <= abs(lit) <= nvariables
               for lit in clause):
            raise ValueError('literal outside declared variable range')
        a, b = node(clause[0]), node(clause[-1])
        for u, v in ((a ^ 1, b), (b ^ 1, a)):
            graph[u].append(v)
            reverse[v].append(u)
    visited = [False] * len(graph)
    finish = []
    for start in range(len(graph)):
        if visited[start]:
            continue
        visited[start] = True
        stack = [(start, 0)]
        while stack:
            vertex, index = stack[-1]
            if index == len(graph[vertex]):
                finish.append(vertex)
                stack.pop()
            else:
                target = graph[vertex][index]
                stack[-1] = vertex, index + 1
                if not visited[target]:
                    visited[target] = True
                    stack.append((target, 0))
    component = [-1] * len(graph)
    label = 0
    for start in reversed(finish):
        if component[start] != -1:
            continue
        component[start] = label
        stack = [start]
        while stack:
            for target in reverse[stack.pop()]:
                if component[target] == -1:
                    component[target] = label
                    stack.append(target)
        label += 1
    if any(component[2 * v] == component[2 * v + 1]
           for v in range(nvariables)):
        return None
    return tuple(int(component[2 * v] > component[2 * v + 1])
                 for v in range(nvariables))


def find_witness(peptide, edges, fixed=None):
    """Return a compatible synonymous RNA or None, in O(m+|H|) RAM operations."""
    variables, clauses = compile_instance(peptide, edges, fixed)
    bits = solve_cnf(variables, clauses)
    if bits is None:
        return None
    return ''.join(BASE[bits[p:p + 2]] for p in range(0, len(bits), 2))
