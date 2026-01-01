"""Independent exact integer characteristic polynomials of feasibility matrices."""
from pathlib import Path
import json
import sympy as sp
ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'exact_peptide_automata.json').read_text());x=sp.Symbol('x');out={}
for d,r in data.items():
    keys=list(r['transition_matrix']);matrix=sp.Matrix([[r['transition_matrix'][i].get(j,0) for j in keys] for i in keys])
    polynomial=matrix.charpoly(x).as_poly()
    out[d]={'degree':polynomial.degree(),'integer_coefficients_descending':[str(c) for c in polynomial.all_coeffs()],
        'factorization':str(sp.factor(polynomial.as_expr()))}
(ROOT/'characteristic_polynomials.json').write_text(json.dumps(out,indent=2))
print('integer characteristic polynomials saved')
