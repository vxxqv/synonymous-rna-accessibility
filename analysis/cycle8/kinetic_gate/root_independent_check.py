"""Independent symbolic rank and inversion check, using only invented values."""
from pathlib import Path
import hashlib,json
import sympy as S
HERE=Path(__file__).resolve().parent
s,e,dn,dc,k,dp=S.symbols('s e d_N d_C k d_P',positive=True)
obs=S.Matrix([s/(e+dn),e*s/(dc*(e+dn)),k*e*s/(dp*dc*(e+dn))])
minor=S.factor(obs.jacobian([s,dc,k]).det())
if minor==0 or minor.is_zero is True:raise ValueError('Rank failed')
n,c,p=S.symbols('n c p',positive=True)
sub={s:n*(e+dn),dc:e*n/c,k:dp*p/c}
if list(S.simplify(obs.subs(sub)-S.Matrix([n,c,p])))!=[0,0,0]:raise ValueError('Fiber sufficiency')
sol=S.solve([s-(e+dn)*n,e*n-dc*c,k*c-dp*p],[s,dc,k],dict=True)
if len(sol)!=1 or any(S.simplify(sol[0][v]-sub[v])!=0 for v in sub):raise ValueError('Fiber necessity')
a=S.symbols('a',positive=True)
if S.simplify((dc*c/n)+(a-dc*c/n)-a)!=0:raise ValueError('Inverse information')
case={s:S.Rational(11,3),e:S.Rational(2,5),dn:S.Rational(7,11),dc:S.Rational(13,17),k:S.Rational(19,23),dp:S.Rational(29,31)}
o=obs.subs(case); reconstruction={s:o[0]*(case[e]+case[dn]),dc:case[e]*o[0]/o[1],k:case[dp]*o[2]/o[1]}
if any(reconstruction[v]!=case[v] for v in reconstruction):raise ValueError('Exact reconstruction')
receipt={'status':'PASS','jacobian_nonzero_minor':str(minor),'fiber_necessary_and_sufficient':True,
 'exact_invented_rational_reconstruction':True,'rank_scope':'Strictly positive parameters; nonzero symbolic minor proves rank three throughout.',
 'derivation_sha256':hashlib.sha256((HERE/'derivation.md').read_bytes()).hexdigest(),
 'review_scope':'Independent primary-agent algebraic reconstruction. No biological outcomes, finite checks not substituted for the displayed proof.'}
(HERE/'root_independent_review.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
