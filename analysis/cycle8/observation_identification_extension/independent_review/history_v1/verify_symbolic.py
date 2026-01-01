"""Independent global symbolic review of fixed equilibrium observations.

No author numerical module is imported. No phenotype or sequence data are read.
"""
import argparse, ast, hashlib, itertools, json, platform, re, time
from fractions import Fraction as Q
from pathlib import Path
import sympy as sp

RATE_SHA='72d2b9d519ac375b55a70ca573501b53bcf440dc6fb9bfbd8df1556aa94af3dc'
PROCESS_SHA='f2a387383c8ed130b1b4f6a846b07c9364bc8eb182668ea578c3165339145ee2'
CHECKS=0
def require(value,reason):
    global CHECKS
    if not value:raise AssertionError(reason)
    CHECKS+=1
def same(u,v):return sp.cancel(u-v)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n',encoding='utf-8')

def global_rate_review(author_inventory):
    # Independent coordinates: total nuclear exit a, fraction f, protein loss z.
    # Domain is a,z,n,c,p>0 and 0<f<1, imposed in the proof rather than inferred
    # from a numerical Jacobian or finite anchor checks.
    a,f,z,n,c,p,tau=sp.symbols('a f z n c p tau',positive=True)
    variables=(a,f,z)
    rates={'s':n*a,'e':a*f,'d_N':a*(1-f),'d_C':n*a*f/c,'k':p*z/c,'d_P':z}
    obs={**rates,'a':a,'J':n*a*f,'H':p*z,'f':f}
    labels=sorted(obs)
    require(len(labels)==10,'Observation list')
    # Equilibrium identities hold symbolically, not only at example anchors.
    require(same(rates['s']-(rates['e']+rates['d_N'])*n,0),'N balance')
    require(same(rates['e']*n-rates['d_C']*c,0),'C balance')
    require(same(rates['k']*c-rates['d_P']*p,0),'P balance')
    inverses=[];zero_determinants=[]
    for triple in itertools.combinations(labels,3):
        jac=sp.Matrix([obs[x] for x in triple]).jacobian(variables);det=sp.factor(jac.det())
        if det==0:zero_determinants.append(triple);continue
        u=sp.symbols('u0:3',positive=True)
        sol=sp.solve([obs[k]-v for k,v in zip(triple,u)],variables,dict=True)
        require(len(sol)==1 and set(sol[0])==set(variables),'Unexpected inverse branches')
        inv=sol[0]
        # Both compositions are checked as rational identities. The report gives
        # the exact positive image restrictions that make each inverse admissible.
        for var in variables:
            require(same(inv[var].subs(dict(zip(u,[obs[k] for k in triple])),simultaneous=True),var),'Inverse composed with observation')
        for key,value in zip(triple,u):require(same(obs[key].subs(inv,simultaneous=True),value),'Observation composed with inverse')
        inverses.append({'observations':triple,'jacobian_determinant':str(det),'inverse':{str(k):str(v) for k,v in inv.items()},
          'inverse_denominators':{str(k):str(sp.factor(sp.denom(sp.cancel(v)))) for k,v in inv.items()}})
    # These maps are globally positive and nontrivial on the declared open domain.
    # Their invariant masks are computed from the actual symbolic observations.
    moves={
      'change_protein_clock':{a:a,f:f,z:2*z},
      'scale_nuclear_total':{a:2*a,f:f,z:z},
      'fixed_total_change_fraction':{a:a,f:f/2,z:z},
      'fixed_export_change_total':{a:2*a,f:f/2,z:z},
      'fixed_nuclear_loss_change_total':{a:2*a,f:(1+f)/2,z:z}}
    invariants={};move_checks={}
    for name,m in moves.items():
        invariants[name]=[label for label in labels if same(obs[label].subs(m,simultaneous=True),obs[label])]
        require(any(not same(var,m[var]) for var in variables),'Nontrivial counterexample family')
        move_checks[name]={'coordinate_map':{str(k):str(v) for k,v in m.items()},'unchanged_observations':invariants[name],
          'positive_domain_proof':'a and z are multiplied by positive constants; f is unchanged, f/2, or (1+f)/2, all strictly inside (0,1).'}
    triples=[set(r['observations']) for r in inverses];classification=[]
    for bits in range(1<<10):
        chosen={labels[i] for i in range(10) if bits&(1<<i)}
        witness=next((i for i,s in enumerate(triples) if s<=chosen),None)
        counter=next((k for k,s in invariants.items() if chosen<=set(s)),None)
        require((witness is None)!=(counter is None),'Every subset needs exactly one global proof type')
        classification.append({'observations':sorted(chosen),'globally_sufficient':witness is not None,
          'minimal':witness is not None and len(chosen)==3,'inverse_certificate_index':witness,'positive_counterexample_family':counter})
    sufficient=sum(x['globally_sufficient'] for x in classification)
    require(len(inverses)==51 and sufficient==805,'Derived classification totals')
    old={frozenset(r['observations']):r for r in author_inventory['subsets']}
    require(len(old)==1024,'Author subset coverage')
    for row in classification:
        saved=old[frozenset(row['observations'])]
        require((saved['globally_sufficient'],saved['minimal'])==(row['globally_sufficient'],row['minimal']),'Author classification disagreement')
    require({frozenset(r['observations']) for r in author_inventory['minimal_sets']}=={frozenset(x) for x in triples},'Minimal triple inventory disagreement')
    # Unknown residence O=p*z*tau: every original ambiguity can preserve occupancy.
    O=p*z*tau;occupancy_maps={}
    for name,m in moves.items():
        aug={**m,tau:tau/2 if name=='change_protein_clock' else tau}
        require(same(O.subs(aug,simultaneous=True),O),'Unknown residence did not preserve occupancy')
        occupancy_maps[name]={str(k):str(v) for k,v in aug.items()}
    require(same(O.subs(tau,O/(p*z)),O),'Residence inverse')
    occupancy_counts={'rate_identifying_among_2048_subsets':0,'augmented_rate_and_residence_identifying_among_2048_subsets':0}
    for row in classification:
        for hasO in [False,True]:
            occupancy_counts['rate_identifying_among_2048_subsets']+=row['globally_sufficient']
            occupancy_counts['augmented_rate_and_residence_identifying_among_2048_subsets']+=row['globally_sufficient'] and hasO
    require(occupancy_counts=={'rate_identifying_among_2048_subsets':1610,'augmented_rate_and_residence_identifying_among_2048_subsets':805},'Occupancy enumeration')
    # Exact derivative separates output-log sensitivity from relative-input elasticity.
    D,F=sp.symbols('D F',positive=True);e=D*F/(1-F)
    dlog=sp.cancel(sp.diff(e,F)/e);elasticity=sp.cancel(F*dlog)
    require(same(dlog,1/(F*(1-F))) and same(elasticity,1/(1-F)),'Conditioning distinction')
    # A real boundary counterexample: D=0,F=1,B=1 does not fix positive export.
    boundary=[{'e':Q(1),'d_N':Q(0),'d_P':Q(1)},{'e':Q(2),'d_N':Q(0),'d_P':Q(1)}]
    require(boundary[0]!=boundary[1] and all(x['e']/(x['e']+x['d_N'])==1 for x in boundary),'Boundary family')
    # Calibrated amounts are necessary, even with an otherwise sufficient triple.
    calibration=[{'n':Q(1),'c':Q(1),'p':Q(1),'n_gain':Q(1),'s':Q(4),'e':Q(1),'d_N':Q(3),'d_C':Q(1),'k':Q(1),'d_P':Q(1)},
                 {'n':Q(2),'c':Q(1),'p':Q(1),'n_gain':Q(1,2),'s':Q(4),'e':Q(1),'d_N':Q(1),'d_C':Q(2),'k':Q(1),'d_P':Q(1)}]
    for t in calibration:
        require(t['s']==(t['e']+t['d_N'])*t['n'] and t['e']*t['n']==t['d_C']*t['c'] and t['k']*t['c']==t['d_P']*t['p'],'Calibration example balance')
    require(calibration[0]['n_gain']*calibration[0]['n']==calibration[1]['n_gain']*calibration[1]['n']==1,'Equal observed nuclear amount')
    require(all(calibration[0][k]==calibration[1][k] for k in ['s','e','d_P','c','p']),'Equal sufficient observations when calibration is missing')
    return {'coordinates':['a=e+d_N','f=e/(e+d_N)','z=d_P'],'domain':'a,z,n,c,p>0; 0<f<1',
      'symbolic_inverse_triples':inverses,'singular_triples':zero_determinants,'counterexample_families':move_checks,
      'classification':classification,'sufficient_subsets':sufficient,'insufficient_subsets':1024-sufficient,
      'occupancy_preserving_counterexample_maps':occupancy_maps,'occupancy_counts':occupancy_counts,
      'conditioning':{'d_log_e_d_f':str(dlog),'d_log_e_d_log_f':str(elasticity)},
      'excluded_boundary_counterexample':boundary,'unknown_calibration_counterexample':calibration}

def processing_review(parser_path,process_path,author_examples):
    A=sp.Matrix([[1,0,0,0],[0,1,1,0]]);ell=sp.Matrix([[1,1,0,0]])
    require(A.rank()==2,'A rank');null=A.nullspace();require(len(null)==2,'A nullity')
    require(A.col_join(ell).rank()==3,'Coding contrast incorrectly row estimable')
    total=A.col_join(sp.Matrix([[1,1,1,1]]))
    require(total.rank()==3 and total.col_join(ell).rank()==4,'Total resolves ambiguity unexpectedly')
    # Independent positive pair differs from the author's example.
    u=sp.Matrix([3,1,4,2]);v=sp.Matrix([3,3,2,2])
    require(A*u==A*v==sp.Matrix([3,5]) and sum(u)==sum(v)==10,'Positive ambiguity observation')
    require((ell*u)[0]==4 and (ell*v)[0]==6,'Positive ambiguity coding contrast')
    # Every binary linear contrast is classified independently by exact row rank.
    contrast_classes=[]
    for bits in itertools.product([0,1],repeat=4):
        l=sp.Matrix([bits]);rowspace=A.col_join(l).rank()==A.rank()
        require(rowspace==(bits[1]==bits[2] and bits[3]==0),'Row-space characterization')
        withtotal=total.col_join(l).rank()==total.rank()
        require(withtotal==(bits[1]==bits[2]),'Row-space with total')
        if not rowspace:
            w=next(w for w in null if (l*w)[0]!=0);x=sp.ones(4,1);xp=x+w/2;xm=x-w/2
            require(all(z>0 for z in list(xp)+list(xm)) and A*xp==A*xm and l*xp!=l*xm,'Kernel positivity counterexample')
        contrast_classes.append({'contrast':bits,'determined_by_A':rowspace,'determined_by_A_and_total':withtotal})
    # Inspect and execute ONLY the original parser function on symbolic predicates.
    tree=ast.parse(parser_path.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='categorize_read')
    required_names={'donor','donor_loc','acceptor','acceptor_loc'}
    assigned={n.id for n in ast.walk(fn) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
    require(required_names<=assigned,'Source parser structure')
    ns={'re':re};exec(compile(ast.Module(body=[fn],type_ignores=[]),'isolated_source_parser','exec'),ns)
    class ReferencePredicates:
        def __getitem__(self,sl):
            t=(sl.start,sl.stop)
            if t in [(10,12),(100,102)]:return 'GT'
            if t in [(38,40),(128,130)]:return 'AG'
            if t==(94,100):return 'terminal_donor'
            if t==(130,136):return 'terminal_acceptor'
            return 'symbolic_other'
    outputs=[ns['categorize_read'](c,d,ReferencePredicates()) for c,d in [('100M30D50M',30),('10M30D60M30D50M',60)]]
    require(outputs[0]==outputs[1]=='AG_junction_low_mismatch-terminal_donor_100_terminal_acceptor_130','Parser compression claim')
    # Source-only AST evidence for the operator. No notebook outputs are read.
    code=process_path.read_text();pt=ast.parse(code)
    functions={n.name:n for n in pt.body if isinstance(n,ast.FunctionDef)}
    require(set(['process_raw_count','filter_barcode_count','compute_median_splice_fraction'])<=set(functions),'Processing source functions')
    med=functions['compute_median_splice_fraction'];calls=[n for n in ast.walk(med) if isinstance(n,ast.Call)]
    pivot=next(n for n in calls if isinstance(n.func,ast.Attribute) and n.func.attr=='pivot_table')
    require('aggfunc' not in {x.arg for x in pivot.keywords},'Unexpected explicit pivot aggregation')
    require(any(isinstance(n.func,ast.Attribute) and n.func.attr=='median' for n in calls),'Median absent')
    filt=functions['filter_barcode_count'];ftext=ast.unparse(filt)
    require('value_counts' in ftext and 'nunique' not in ftext,'Distinct-barcode claim changed')
    # Exact original and independent alternative normalized-median counterexamples.
    def med(values):
        x=sorted(values);return (x[2]+x[3])/2
    counts1=[(80,10,10)]*2+[(10,80,10)]*2+[(10,10,80)]*2
    counts2=[(60,30,10)]*2+[(10,60,30)]*2+[(30,10,60)]*2
    medians=[]
    for counts in [counts1,counts2]:
        require(all(sum(row)==100 and row[0]>0 for row in counts),'Synthetic eligible count design')
        ms=[med([Q(row[j],sum(row)) for row in counts]) for j in range(3)]
        medians.append(ms)
    require(medians[0]==[Q(1,10)]*3 and sum(medians[0])==Q(3,10),'Original median counterexample')
    require(medians[1]==[Q(3,10)]*3 and sum(medians[1])==Q(9,10),'Independent median counterexample')
    require(3*2>5 and 3<=5,'Concatenated full-length rows vs distinct barcodes')
    present_S=Q(1,10);zero_completed_S=(Q(1,10)+0)/2;Fmean=(Q(9,10)+1)/2
    require(present_S!=zero_completed_S and Fmean+present_S==Q(21,20),'Missing-row mean counterexample')
    old_u=sp.Matrix([sp.Rational(x) for x in author_examples['mixtures'][0]]);old_v=sp.Matrix([sp.Rational(x) for x in author_examples['mixtures'][1]])
    require(A*old_u==A*old_v==sp.Matrix([2,4]) and (ell*old_u)[0]==5 and (ell*old_v)[0]==3,'Author positive mixture example')
    return {'known_A_rank':2,'known_A_nullity':2,'coding_contrast_in_rowspace':False,'coding_contrast_in_rowspace_with_total':False,
      'independent_positive_mixtures':[list(u),list(v)],'same_observation':[3,5],'same_total':10,'coding_values':[4,6],
      'binary_contrast_classes':contrast_classes,'isolated_parser_outputs':outputs,
      'source_function_line_ranges':{'categorize_read':[fn.lineno,fn.end_lineno],**{k:[v.lineno,v.end_lineno] for k,v in functions.items()}},
      'normalized_barcode_counterexample_medians':medians,'missing_preparation_example':{'full_length_mean':Fmean,'present_only_spliced_mean':present_S,'zero_completed_spliced_mean':zero_completed_S},
      'no_actual_category_counts_or_read_sequences_accessed':True}

def main():
    if not __debug__:raise RuntimeError('Optimized execution prohibited')
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,required=True);ap.add_argument('--dao-gate',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    paths={'rate':a.source_dir/'rate_observation_theorem.md','processing':a.source_dir/'processing_observation_contract.md',
      'subsets':a.source_dir/'observation_subsets.json','examples':a.source_dir/'rational_examples.json',
      'parser':a.dao_gate/'source_junction_parser.py','operator':a.dao_gate/'processing_code_only.txt'}
    pins={k:sha(p) for k,p in paths.items()}
    require(pins['rate']==RATE_SHA and pins['processing']==PROCESS_SHA,'Held claim hash mismatch')
    started=time.perf_counter();a.out.mkdir(parents=True)
    try:
        G=global_rate_review(json.loads(paths['subsets'].read_text()))
        C=processing_review(paths['parser'],paths['operator'],json.loads(paths['examples'].read_text()))
        require(pins=={k:sha(p) for k,p in paths.items()},'Source drift')
        put(a.out/'global_symbolic_evidence.json',G);put(a.out/'processing_evidence.json',C)
        receipt={'status':'PASS_WITH_CONDITIONING_WORDING_CLARIFICATION','checks':CHECKS,'rate_subsets':1024,'global_inverse_triples':len(G['symbolic_inverse_triples']),
          'sufficient':G['sufficient_subsets'],'insufficient':G['insufficient_subsets'],'input_sha256':pins,
          'source_sha256':sha(Path(__file__)),'output_sha256':{n:sha(a.out/n) for n in ['global_symbolic_evidence.json','processing_evidence.json']},
          'runtime':{'python':platform.python_version(),'sympy':sp.__version__},'elapsed_seconds':time.perf_counter()-started,
          'observed_outcomes_read':False,'new_nucleotide_sequences_created':False,'limitation':'Global proofs use the explicitly positive calibrated equilibrium model. Source parser test uses symbolic slice predicates only.'}
        put(a.out/'receipt.json',receipt);print(json.dumps(receipt))
    except Exception as exc:
        put(a.out/'failure.json',{'status':'BLOCKED','error':str(exc),'checks':CHECKS,'input_sha256':pins});raise

if __name__=='__main__':main()
