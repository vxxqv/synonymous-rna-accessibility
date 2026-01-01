"""Independent B1 reference. No production module is imported."""
from fractions import Fraction as Q
from collections import defaultdict, deque
from pathlib import Path
from copy import deepcopy
from io import BytesIO
import hashlib,json,zipfile,xml.etree.ElementTree as E

NS='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
N={'m':NS}
REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
CELL=['K562','KBM7']
GAMMAS=[Q(1),Q(5,4),Q(3,2),Q(2)]
GRID=[Q(5,4),Q(3,2),Q(2),Q(4)]

def fraction(value):
    if value is None or isinstance(value,bool):return None
    try:return Q(str(value))
    except (ValueError,TypeError,ZeroDivisionError):return None

def ready(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {str(k):ready(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [ready(v) for v in x]
    return x

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def read_selected(source,eligibility):
    """XML-select first, then openpyxl-read a sanitized in-memory workbook.

    Unselected numeric cells are never parsed as numbers. Exact source decimals
    supply arithmetic; openpyxl is a separate paired-cell/roundtrip check.
    """
    import openpyxl
    z=zipfile.ZipFile(source);wb=E.fromstring(z.read('xl/workbook.xml'))
    ids={s.attrib['name']:s.attrib['{'+REL+'}id'] for s in wb.findall('m:sheets/m:sheet',N)}
    links=E.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    target=next(r.attrib['Target'] for r in links if r.attrib['Id']==ids['Supplementary_Table_3_Multi_exo'])
    target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
    ss=E.fromstring(z.read('xl/sharedStrings.xml')) if 'xl/sharedStrings.xml' in z.namelist() else None
    cache={}
    def text(c):
        if c is None:return None
        kind=c.attrib.get('t');v=c.find('m:v',N)
        if kind=='inlineStr':return ''.join(a.text or '' for a in c.findall('.//m:t',N))
        if v is None:return None
        if kind=='s':
            i=int(v.text)
            if i not in cache:cache[i]=''.join(a.text or '' for a in ss[i].findall('.//m:t',N))
            return cache[i]
        if kind=='b':return '[BOOLEAN]'
        if kind=='e':return '[EXCEL_ERROR]'
        return v.text
    tree=E.fromstring(z.read(target));sheet=tree.find('m:sheetData',N)
    allowed=eligibility['states'];selected={};safe=deepcopy(tree);safe_data=safe.find('m:sheetData',N)
    for child in list(safe_data):safe_data.remove(child)
    for row in sheet:
        rownum=int(row.attrib['r']);cells={''.join(x for x in c.attrib['r'] if x.isalpha()):c for c in row}
        if rownum==1:
            expected={'A':'combined_var_id','B':'var_id','E':'exon_number','F':'syn_var','G':'cell_line','H':'RPM_DNA','I':'RPM_RNA','J':'cprime_score'}
            if any(text(cells.get(c))!=v for c,v in expected.items()):raise ValueError('Changed header schema')
            clean=deepcopy(row)
            for c in list(clean):
                if ''.join(x for x in c.attrib['r'] if x.isalpha()) not in list('ABCDEFGHIJ'):clean.remove(c)
            safe_data.append(clean);continue
        state=text(cells.get('A'))
        if state not in allowed:continue
        cell=text(cells.get('G'))
        if cell not in eligibility['cell_contexts']:raise ValueError('Unexpected selected context')
        meta=allowed[state]
        if text(cells.get('B'))!=meta['var_id'] or text(cells.get('F'))!=meta['syn_var']:raise ValueError('Metadata identity mismatch')
        exon=fraction(text(cells.get('E')))
        if exon!=Q(str(meta['exon_number'])):raise ValueError('Exon mismatch')
        values={name:text(cells.get(col)) for name,col in [('DNA','H'),('RNA','I'),('source_score','J')]}
        selected[rownum]={'source_row':rownum,'state':state,'cell':cell,'intended':meta['var_id'],'marker':meta['syn_var'],'exon':int(exon),'raw':values}
        clean=deepcopy(row)
        for c in list(clean):
            col=''.join(x for x in c.attrib['r'] if x.isalpha())
            if col not in list('ABCDEFGHIJ'):clean.remove(c)
            elif col in ['H','I','J'] and c.attrib.get('t','n')=='n' and fraction(text(c)) is None:
                # Invalid numeric XML becomes a harmless string only in the
                # openpyxl cross-check copy; the original missingness is retained.
                original=text(c);c.attrib['t']='inlineStr'
                for child in list(c):c.remove(child)
                E.SubElement(E.SubElement(c,'{'+NS+'}is'),'{'+NS+'}t').text=original or ''
        safe_data.append(clean)
    counter=defaultdict(int)
    for r in selected.values():counter[r['state'],r['cell']]+=1
    expected={(s,c) for s in allowed for c in eligibility['cell_contexts']}
    if set(counter)!=expected or any(n!=3 for n in counter.values()):raise ValueError('Selected metadata multiplicity mismatch')
    # openpyxl sees no unselected response cell, including poisoned synthetic cells.
    buf=BytesIO()
    with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as out:
        for name in z.namelist():out.writestr(name,E.tostring(safe,encoding='utf-8',xml_declaration=True) if name==target else z.read(name))
    buf.seek(0);w=openpyxl.load_workbook(buf,read_only=True,data_only=True)
    s=w['Supplementary_Table_3_Multi_exo'];roundtrip=[];rows=[]
    for n,row in enumerate(s.iter_rows(min_col=8,max_col=10,values_only=True),1):
        if n not in selected:continue
        r=selected[n]
        for k,v in zip(['DNA','RNA','source_score'],row):
            exact=fraction(r['raw'][k]);rounded=fraction(v)
            if exact!=rounded:roundtrip.append({'source_row':n,'field':k,'source_decimal':r['raw'][k],'openpyxl_repr':repr(v),'exact_difference':None if exact is None or rounded is None else exact-rounded})
        rows.append(row_result(r))
    return rows,roundtrip

def row_result(r):
    vals={k:fraction(v) for k,v in r['raw'].items()}
    dna,rna,score=[vals[k] for k in ['DNA','RNA','source_score']]
    ratio=rna/dna if dna is not None and rna is not None and dna>0 and rna>0 else None
    delta=score-ratio if ratio is not None and score is not None else None
    consistent=abs(delta)<=Q(1,10**8)*max(1,abs(ratio)) if delta is not None else None
    reasons=[k for k,v in vals.items() if v is None or v<=0]
    return {**{k:v for k,v in r.items() if k!='raw'},**vals,'ratio':ratio,
            'score_difference':delta,'source_score_consistent':consistent,
            'usable':not reasons,'missing_reason':reasons}

def states_from_rows(rows):
    by=defaultdict(list)
    for r in rows:by[r['state'],r['cell']].append(r)
    out=[]
    for (state,cell),rr in sorted(by.items()):
        if len(rr)!=3:raise ValueError('Each state/cell requires three rows')
        values=sorted(r['ratio'] for r in rr if r['usable'])
        complete=len(values)==3
        out.append({'state':state,'cell':cell,'usable_rows':len(values),'complete':complete,
                    'lower':values[0] if complete else Q(0),'upper':values[-1] if complete else None,
                    'median':values[1] if complete else None,'source_rows':sorted(r['source_row'] for r in rr)})
    return out

def group_bound(intended,cell,markers,states,gamma):
    if len(markers)<2:raise ValueError('Marker contrast requires two states')
    boxes=[(states[s,cell]['lower']/gamma,None if states[s,cell]['upper'] is None else states[s,cell]['upper']*gamma) for s in markers]
    finite_upper=[b for a,b in boxes if b is not None]
    lo=max(a for a,b in boxes);hi=min(finite_upper) if finite_upper else None
    minimum=max(Q(1),lo/hi) if hi is not None else Q(1)
    infinite=any(b is None or a==0 for a,b in boxes)
    maximum=None if infinite else max([Q(1)]+[boxes[i][1]/boxes[j][0] for i in range(len(boxes)) for j in range(len(boxes)) if i!=j])
    centers=[states[s,cell]['median'] for s in markers]
    point=None if None in centers else max(centers)/min(centers)
    return {'intended':intended,'cell':cell,'markers':sorted(markers),'marker_count':len(markers),
            'complete':all(states[s,cell]['complete'] for s in markers),
            'minimum_fold':minimum,'maximum_fold':maximum,'maximum_attained':maximum is not None,
            'point_fold':point,'common_response_feasible':hi is None or lo<=hi,
            'thresholds':{str(t):{'minimum_at_least':minimum>=t,'maximum_at_most':maximum is not None and maximum<=t} for t in GRID}}

def shared_order(exon,eligible,states,gamma):
    nodes=sorted(s for s,m in eligible['states'].items() if int(m['exon_number'])==exon)
    boxes={(s,c):(states[s,c]['lower']/gamma,None if states[s,c]['upper'] is None else states[s,c]['upper']*gamma) for s in nodes for c in eligible['cell_contexts']}
    edges={};adj={s:set() for s in nodes};indegree=dict.fromkeys(nodes,0)
    for a in nodes:
        for b in nodes:
            if a==b:continue
            origins=[{'context':c,'upper':boxes[a,c][1],'lower':boxes[b,c][0]} for c in eligible['cell_contexts'] if boxes[a,c][1] is not None and boxes[a,c][1]<boxes[b,c][0]]
            if origins:edges[a,b]={'lower_state':a,'higher_state':b,'origins':origins};adj[a].add(b);indegree[b]+=1
    queue=deque(sorted(s for s in nodes if indegree[s]==0));order=[]
    while queue:
        a=queue.popleft();order.append(a)
        for b in sorted(adj[a]):
            indegree[b]-=1
            if indegree[b]==0:queue.append(b)
    out={'exon':exon,'state_count':len(nodes),'status':'feasible' if len(order)==len(nodes) else 'infeasible','edges':[edges[k] for k in sorted(edges)]}
    if len(order)==len(nodes):
        responses={}
        for c in eligible['cell_contexts']:
            upper=[u for s in nodes for l,u in [boxes[s,c]] if u is not None]
            last=min([Q(1)]+upper)/2;values={}
            for s in order:
                l,u=boxes[s,c];last=max(last,l)
                if u is not None and last>u:raise ArithmeticError('Invalid independent feasible completion')
                values[s]=last
            responses[c]=values
        out.update(latent_order=order,responses=responses)
    else:
        color={s:0 for s in nodes};stack=[]
        def visit(a):
            color[a]=1;stack.append(a)
            for b in sorted(adj[a]):
                if color[b]==1:return stack[stack.index(b):]+[b]
                if color[b]==0:
                    found=visit(b)
                    if found:return found
            stack.pop();color[a]=2
        cycle=None
        for a in nodes:
            if color[a]==0:
                cycle=visit(a)
                if cycle:break
        if not cycle:raise ArithmeticError('No cycle found in cyclic graph')
        out.update(cycle=cycle,cycle_edges=[edges[a,b] for a,b in zip(cycle,cycle[1:])])
    return out

def reconstruct(rows,eligible):
    bad=[r['source_row'] for r in rows if r['source_score_consistent'] is False]
    if bad:return {'status':'SOURCE_REPRODUCTION_FAILED','source_rows':rows,'inconsistent_rows':bad}
    state_list=states_from_rows(rows);states={(s['state'],s['cell']):s for s in state_list};sens=[]
    for gamma in GAMMAS:
        groups=[group_bound(i,c,markers,states,gamma) for i,markers in sorted(eligible['groups'].items()) for c in eligible['cell_contexts']]
        orders=[shared_order(exon,eligible,states,gamma) for exon in eligible['exons']]
        summary={'group_cells':len(groups),'complete_group_cells':sum(g['complete'] for g in groups),
                 'strictly_positive_lower_dispersion':sum(g['minimum_fold']>1 for g in groups),
                 'minimum_at_least_twofold':sum(g['minimum_fold']>=2 for g in groups),
                 'maximum_at_most_twofold':sum(g['maximum_fold'] is not None and g['maximum_fold']<=2 for g in groups),
                 'shared_order_infeasible_exons':sum(o['status']=='infeasible' for o in orders)}
        sens.append({'gamma':gamma,'group_cell_results':groups,'exon_shared_order_results':orders,'summary':summary})
    return {'status':'RECONSTRUCTED','source_rows':rows,'state_cells':state_list,'sensitivities':sens}
