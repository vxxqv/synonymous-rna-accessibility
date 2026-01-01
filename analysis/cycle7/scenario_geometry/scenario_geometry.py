"""Exact finite-scenario path quantifiers. No file reader or observed-data mode."""
from fractions import Fraction
from itertools import combinations
import math

LABELS=('fixed_route','adaptive_routes_only','mixed_scenario_connectivity',
        'joint_obstruction_box_connected','upper_box_obstruction')
METRICS=('E','C_fix','C_ad','C_max','C_upper')
NEG_INF=-math.inf
POS_INF=math.inf

def require(condition,message):
    if not condition:raise ValueError(message)

def value(x):
    if x is None:return None
    require(type(x) is not bool and isinstance(x,(int,str,Fraction)), 'Defined scores must be exact rationals')
    q=Fraction(x)
    require(q>=0,'Scores must be nonnegative')
    return q

def normalized(fields,edges):
    rows=[tuple(value(x) for x in row) for row in fields]
    require(bool(rows) and bool(rows[0]),'At least one nonempty field is required')
    n=len(rows[0])
    require(all(len(row)==n for row in rows),'Ragged fields')
    pairs=[];seen=set()
    for edge in edges:
        require(len(edge)==2,'Edge arity')
        a,b=edge
        require(type(a) is int and type(b) is int and 0<=a<n and 0<=b<n and a!=b,'Invalid edge')
        pair=tuple(sorted((a,b)))
        require(pair not in seen,'Duplicate edge')
        seen.add(pair);pairs.append(pair)
    return rows,tuple(sorted(pairs))

def widest(scores,edges):
    """Max-min Floyd closure, with endpoint scores and disconnected value -inf."""
    n=len(scores)
    a=[[NEG_INF]*n for _ in range(n)]
    for i,x in enumerate(scores):a[i][i]=x
    for i,j in edges:a[i][j]=a[j][i]=min(scores[i],scores[j])
    for k in range(n):
        for i in range(n):
            if a[i][k]==NEG_INF:continue
            for j in range(n):
                candidate=min(a[i][k],a[k][j])
                if candidate>a[i][j]:a[i][j]=candidate
    return a

def capacity_tables(fields,edges):
    n=len(fields[0])
    per_scenario=[widest(row,edges) for row in fields]
    low=[min(row[i] for row in fields) for i in range(n)]
    high=[max(row[i] for row in fields) for i in range(n)]
    fixed=widest(low,edges);upper=widest(high,edges)
    records=[]
    for s,t in combinations(range(n),2):
        caps=tuple(matrix[s][t] for matrix in per_scenario)
        record={'s':s,'t':t,'scenario_capacities':caps,
                'E':min(low[s],low[t]),'C_fix':fixed[s][t],
                'C_ad':min(caps),'C_max':max(caps),'C_upper':upper[s][t]}
        require(record['C_fix']<=record['C_ad']<=record['C_max']<=record['C_upper'],'Capacity order violation')
        require(record['C_ad']<=record['E'],'Endpoint bound violation')
        records.append(record)
    return records

def intersect_interval(label,lower,upper,E):
    """Intersect (lower,upper] with [0,E]; upper infinity is harmless here."""
    lo=max(Fraction(0),lower);hi=min(E,upper)
    closed=lower<0
    nonempty=lo<hi or lo==hi and closed
    return {'label':label,'lower':lo,'upper':hi,'lower_closed':closed,
            'upper_closed':True,'nonempty':nonempty,
            'width':max(Fraction(0),hi-lo)}

def intervals(record):
    cuts=(NEG_INF,record['C_fix'],record['C_ad'],record['C_max'],record['C_upper'],POS_INF)
    return [intersect_interval(label,cuts[i],cuts[i+1],record['E']) for i,label in enumerate(LABELS)]

def classify_exact(record,threshold):
    tau=value(threshold);require(tau is not None,'Threshold cannot be unknown')
    if tau>record['E']:return 'ineligible_endpoints'
    for label,key in zip(LABELS[:-1],('C_fix','C_ad','C_max','C_upper')):
        if tau<=record[key]:return label
    return LABELS[-1]

def analyze(fields,edges):
    fields,edges=normalized(fields,edges)
    unknown=[(r,v) for r,row in enumerate(fields) for v,x in enumerate(row) if x is None]
    if not unknown:
        records=capacity_tables(fields,edges)
        for record in records:
            record['status']='exact'
            record['intervals']=intervals(record)
            record['all_scenario_obstruction']=intersect_interval('all_scenario_obstruction',record['C_max'],POS_INF,record['E'])
            record['upper_box_obstruction']=intersect_interval('upper_box_obstruction',record['C_upper'],POS_INF,record['E'])
    else:
        lower=capacity_tables([[Fraction(0) if x is None else x for x in row] for row in fields],edges)
        upper=capacity_tables([[POS_INF if x is None else x for x in row] for row in fields],edges)
        records=[]
        for low,high in zip(lower,upper):
            require((low['s'],low['t'])==(high['s'],high['t']),'Bound indexing')
            bounds={key:(low[key],high[key]) for key in METRICS}
            require(all(a<=b for a,b in bounds.values()),'Bound order')
            records.append({'s':low['s'],'t':low['t'],'status':'bounded_unknown',
                            'bounds':bounds,'scenario_capacity_bounds':tuple(zip(low['scenario_capacities'],high['scenario_capacities'])),
                            'intervals':None})
    return {'vertices':len(fields[0]),'scenarios':len(fields),'edges':len(edges),
            'unknown_cells':unknown,'pair_count':len(records),'pairs':records}

def classify_bounded(record,threshold):
    tau=value(threshold);require(tau is not None,'Threshold cannot be unknown')
    bounds=record['bounds']
    if tau>bounds['E'][1]:return 'ineligible_endpoints'
    if tau>bounds['E'][0]:return 'unresolved_endpoint_eligibility'
    if tau<=bounds['C_fix'][0]:return LABELS[0]
    if bounds['C_fix'][1]<tau<=bounds['C_ad'][0]:return LABELS[1]
    if bounds['C_ad'][1]<tau<=bounds['C_max'][0]:return LABELS[2]
    if bounds['C_max'][1]<tau<=bounds['C_upper'][0]:return LABELS[3]
    if bounds['C_upper'][1]<tau:return LABELS[4]
    return 'unresolved_classification'

def json_ready(obj):
    if isinstance(obj,Fraction):return str(obj)
    if isinstance(obj,float):
        require(math.isinf(obj),'Unexpected floating value')
        return 'inf' if obj>0 else '-inf'
    if isinstance(obj,dict):return {str(k):json_ready(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):return [json_ready(v) for v in obj]
    return obj
