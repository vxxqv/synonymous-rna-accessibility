"""Explicitly post-primary reporting diagnostic; never changes any result flag."""
if not __debug__:raise RuntimeError('Optimized execution is forbidden.')
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
from collections import Counter
import csv,hashlib,json

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'cprime_results/primary_v1/result.json'
PIN='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'
PROXIMITY=F(1,10**12)


def decimal(x):
    with localcontext() as c:
        c.prec=65
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==PIN
    result=json.loads(SOURCE.read_text())
    allrows=[];cycles=[];summaries=[]
    for s in result['sensitivities']:
        records=[]
        for g in s['group_cell_results']:
            for name,field,t in [('minimum_vs_one','minimum_fold',F(1)),('minimum_vs_two','minimum_fold',F(2)),('maximum_vs_two','maximum_fold',F(2))]:
                assert g[field] is not None
                value=F(g[field]);gap=value-t;relative=gap/t
                sign='equal' if gap==0 else 'above' if gap>0 else 'below'
                row={'gamma':s['gamma'],'intended':g['intended'],'cell':g['cell'],'diagnostic':name,
                     'value_exact':str(value),'value_decimal':decimal(value),'threshold':str(t),
                     'signed_margin_exact':str(gap),'signed_margin_decimal':decimal(gap),
                     'relative_margin_exact':str(relative),'relative_margin_decimal':decimal(relative),
                     'side':sign,'within_postprimary_relative_1e_minus_12':abs(relative)<=PROXIMITY,
                     'frozen_decision':value>t if name=='minimum_vs_one' else value>=t if name=='minimum_vs_two' else value<=t}
                allrows.append(row);records.append(row)
        bydiagnostic={}
        for name in ['minimum_vs_one','minimum_vs_two','maximum_vs_two']:
            items=[x for x in records if x['diagnostic']==name]
            near=[x for x in items if x['within_postprimary_relative_1e_minus_12']]
            positives=[x for x in items if x['frozen_decision']]
            bydiagnostic[name]={'frozen_true_count':len(positives),'near_by_side':dict(Counter(x['side'] for x in near)),
                                'frozen_true_and_near':sum(x['frozen_decision'] for x in near),
                                'frozen_false_and_near':sum(not x['frozen_decision'] for x in near),
                                'nonzero_near_count':sum(x['side']!='equal' for x in near)}
        summaries.append({'gamma':s['gamma'],'diagnostics':bydiagnostic})
        for e in s['exon_shared_order_results']:
            for edge in e.get('cycle_edges',[]):
                for o in edge['origins']:
                    u,l=F(o['upper']),F(o['lower'])
                    gap=l-u;relative=gap/max(abs(l),abs(u))
                    assert gap>0
                    cycles.append({'gamma':s['gamma'],'exon':e['exon'],'lower_state':edge['lower_state'],
                                   'higher_state':edge['higher_state'],'context':o['context'],
                                   'upper_exact':str(u),'lower_exact':str(l),
                                   'gap_exact':str(gap),'gap_decimal':decimal(gap),
                                   'relative_gap_exact':str(relative),'relative_gap_decimal':decimal(relative),
                                   'within_postprimary_relative_1e_minus_12':relative<=PROXIMITY})
    write_csv(HERE/'boundary_proximity_all.csv',allrows)
    nonzero=[x for x in allrows if x['within_postprimary_relative_1e_minus_12'] and x['side']!='equal']
    write_csv(HERE/'boundary_proximity_nonzero_near.csv',nonzero)
    write_csv(HERE/'cycle_margins.csv',cycles)
    report={'status':'postprimary_diagnostic_complete','source_sha256':PIN,'reporting_relative_threshold_exact':str(PROXIMITY),
            'range_margin_definition':'(value-threshold)/threshold, threshold1or2',
            'cycle_margin_definition':'(higher_lower-lower_upper)/max(abs(higher_lower),abs(lower_upper))',
            'all_group_comparisons':len(allrows),'summaries':summaries,
            'nonzero_near_records':nonzero,'stored_cycle_inequalities':len(cycles),
            'cycle_near_count':sum(x['within_postprimary_relative_1e_minus_12'] for x in cycles),
            'smallest_cycle_relative_gap_exact':str(min(F(x['relative_gap_exact']) for x in cycles)),
            'smallest_cycle_relative_gap_decimal':decimal(min(F(x['relative_gap_exact']) for x in cycles)),
            'frozen_results_or_thresholds_changed':False,'statistical_or_biological_error_bound':False,
            'source_rounding_procedure_authenticated':False}
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==PIN
    (HERE/'boundary_proximity.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'summaries':summaries,'nonzero_near_records':nonzero,'stored_cycle_inequalities':len(cycles),
                      'cycle_near_count':report['cycle_near_count'],'smallest_cycle_relative_gap_decimal':report['smallest_cycle_relative_gap_decimal']},indent=2))


if __name__=='__main__':main()
