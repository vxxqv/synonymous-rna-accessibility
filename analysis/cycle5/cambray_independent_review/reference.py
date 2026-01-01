"""Independently authored Cambray protocol arithmetic; no primary imports/I/O."""
import sys
sys.dont_write_bytecode=True
import math
from collections import defaultdict

METHODS=('additive','neighbor_mean','three_mean')
SQUARE_EDGES=((0,1),(0,2),(1,3),(2,3))
DIAGONALS=((0,3),(1,2))

def closure_connected(vertices,edges,allowed,s,t):
    """Boolean transitive closure, independently of any queue-based path code."""
    retained=[v for v in vertices if v in allowed]
    if s not in retained or t not in retained:return False
    index={v:i for i,v in enumerate(retained)};n=len(retained)
    matrix=[[i==j for j in range(n)] for i in range(n)]
    for a,b in edges:
        if a in index and b in index:matrix[index[a]][index[b]]=matrix[index[b]][index[a]]=True
    for k in range(n):
        for i in range(n):
            if matrix[i][k]:
                for j in range(n):matrix[i][j]=matrix[i][j] or matrix[k][j]
    return matrix[index[s]][index[t]]

def finite(value):return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)
def pass_value(value,tau,tolerance=0):return finite(value) and value>=tau-tolerance
def classify(values,edges,s,t,thresholds=None,tolerances=(0,0)):
    if any(not finite(values[v][k]) for v in (s,t) for k in range(2)):
        return {'status':'ineligible_endpoint'}
    if any(len(row)!=2 or any(not finite(x) for x in row) for row in values.values()):
        return {'status':'incomplete_point_field'}
    if thresholds is None:
        thresholds=tuple(min(values[s][k],values[t][k]) for k in range(2))
    if any(not pass_value(values[v][k],thresholds[k],tolerances[k]) for v in (s,t) for k in range(2)):
        return {'status':'ineligible_endpoint'}
    allowed=[{v for v,row in values.items() if pass_value(row[k],thresholds[k],tolerances[k])} for k in range(2)]
    individual=[closure_connected(list(values),edges,allowed[k],s,t) for k in range(2)]
    joint=closure_connected(list(values),edges,allowed[0]&allowed[1],s,t)
    status='joint_connected' if joint else 'joint_only' if all(individual) else 'individual_disconnected'
    return {'status':status,'individual':individual,'joint':joint,'thresholds':list(thresholds)}

def bounds(intervals,edges,s,t,thresholds):
    """Unknown traits are None; one known failure blocks optimistic joint membership."""
    lower=[];upper=[]
    for k in range(2):
        lower.append({v for v,row in intervals.items() if row[k] is not None and row[k][0]>=thresholds[k]})
        upper.append({v for v,row in intervals.items() if row[k] is None or row[k][1]>=thresholds[k]})
    if any(v not in lower[k] for v in (s,t) for k in range(2)):
        return {'status':'ineligible_endpoint'}
    pessimistic=[closure_connected(list(intervals),edges,lower[k],s,t) for k in range(2)]
    optimistic=[closure_connected(list(intervals),edges,upper[k],s,t) for k in range(2)]
    joint_pess=closure_connected(list(intervals),edges,lower[0]&lower[1],s,t)
    joint_opt=closure_connected(list(intervals),edges,upper[0]&upper[1],s,t)
    certificate=all(pessimistic) and not joint_opt
    return {'status':'joint_only_certificate' if certificate else 'not_joint_only_certificate',
            'individual_pessimistic':pessimistic,'individual_optimistic':optimistic,
            'joint_pessimistic':joint_pess,'joint_optimistic':joint_opt,
            'joint_only_certificate':certificate}

def preparation_summary(replicates,edges=SQUARE_EDGES,s=0,t=3):
    """Rows contain four protein values then three growth values, never an assay table."""
    complete=all(len(row)==7 and all(finite(x) for x in row) for row in replicates.values())
    if not complete:return {'status':'incomplete_required_measurements'}
    means={v:(math.fsum(row[:4])/4,math.fsum(row[4:])/3) for v,row in replicates.items()}
    flags=[classify({v:(row[i],row[4+j]) for v,row in replicates.items()},edges,s,t)['status']
           for i in range(4) for j in range(3)]
    intervals={v:((min(row[:4]),max(row[:4])),(min(row[4:]),max(row[4:]))) for v,row in replicates.items()}
    tau=[min(intervals[v][k][0] for v in (s,t)) for k in range(2)]
    return {'status':'complete','means':means,'point':classify(means,edges,s,t),
            'preparation_statuses':flags,'all12_joint_only':all(x=='joint_only' for x in flags),
            'range_thresholds':tau,'range_certificate':bounds(intervals,edges,s,t,tau)}

def interaction(values):
    if len(values)!=4 or not all(finite(v) for v in values):return None
    h=(1,-1,-1,1);contrast=math.fsum(sign*v for sign,v in zip(h,values))
    result={'contrast':contrast,'projection':[v-sign*contrast/4 for v,sign in zip(values,h)]}
    result['natural_log']=interaction_log(values) if min(values)>0 else None
    return result
def interaction_log(values):
    y=[math.log(v) for v in values];h=(1,-1,-1,1);c=math.fsum(sign*v for sign,v in zip(h,y))
    return {'contrast':c,'projection':[v-sign*c/4 for v,sign in zip(y,h)]}

def local_predictions(known,target):
    if target not in range(4) or set(known)!=set(range(4))-{target}:
        raise ValueError('Exactly three authorized corner indices required')
    # Key validation is complete before any value conversion or access.
    if any(len(row)!=2 or any(not finite(v) for v in row) for row in known.values()):
        return None
    a,b,opposite=target^1,target^2,target^3
    return {'additive':[math.fsum([known[a][k],known[b][k],-known[opposite][k]]) for k in range(2)],
            'neighbor_mean':[math.fsum([known[a][k],known[b][k]])/2 for k in range(2)],
            'three_mean':[math.fsum(row[k] for row in known.values())/3 for k in range(2)]}

def canonical_cells(values,cap):
    """Disjoint weak-acceptance cells (lo,hi], with representative hi."""
    boundaries=sorted({v for v in values if finite(v) and v<=cap}|{cap})
    return [(None if i==0 else boundaries[i-1],hi,hi) for i,hi in enumerate(boundaries)]
def critical_regions(values,edges,s,t):
    caps=[min(values[s][k],values[t][k]) for k in range(2)]
    axes=[canonical_cells([row[k] for row in values.values()],caps[k]) for k in range(2)]
    return [{'cells':[a[:2],b[:2]],'representative':[a[2],b[2]],
             'classification':classify(values,edges,s,t,(a[2],b[2]))['status']} for a in axes[0] for b in axes[1]]

def grouped_metrics(cases):
    """Complete metric cases have group, truth(2), predictions(3x2), true/pred class."""
    grouped=defaultdict(list)
    for row in cases:grouped[row['group']].append(row)
    metrics={}
    for method in METHODS:
        if not grouped:
            metrics[method]={'MAE':[None,None],'RMSE':[None,None],'misclassification':None,
                             'confusion':{'TP':0,'FN':0,'FP':0,'TN':0},'balanced_accuracy':None};continue
        mae=[];rmse=[]
        for k in range(2):
            group_abs=[];group_sq=[]
            for rows in grouped.values():
                errors=[row['predictions'][method][k]-row['truth'][k] for row in rows]
                group_abs.append(math.fsum(abs(x) for x in errors)/len(errors))
                group_sq.append(math.fsum(x*x for x in errors)/len(errors))
            mae.append(math.fsum(group_abs)/len(group_abs));rmse.append(math.sqrt(math.fsum(group_sq)/len(group_sq)))
        mis=math.fsum(sum(row['predicted_class'][method]!=row['true_class'] for row in rows)/len(rows)
                      for rows in grouped.values())/len(grouped)
        confusion={'TP':0,'FN':0,'FP':0,'TN':0}
        for row in cases:
            true,pred=row['true_class'],row['predicted_class'][method]
            confusion['TP' if true and pred else 'FN' if true else 'FP' if pred else 'TN']+=1
        p,n=confusion['TP']+confusion['FN'],confusion['TN']+confusion['FP']
        ba=.5*(confusion['TP']/p+confusion['TN']/n) if p and n else None
        metrics[method]={'MAE':mae,'RMSE':rmse,'misclassification':mis,'confusion':confusion,'balanced_accuracy':ba}
    both_classes={row['true_class'] for row in cases}=={False,True}
    coverage=len(cases)>=10 and len(grouped)>=8
    superior=coverage and all(metrics[baseline]['MAE'][k]>0 and
           metrics['additive']['MAE'][k]<=.9*metrics[baseline]['MAE'][k]
           for baseline in METHODS[1:] for k in range(2))
    no_worse=bool(cases) and all(metrics['additive']['misclassification']<=metrics[b]['misclassification'] for b in METHODS[1:])
    return {'n_cases':len(cases),'n_groups':len(grouped),'metrics':metrics,
            'both_classes':both_classes,'coverage_pass':coverage,
            'numeric_superiority':superior,'classification_no_worse':no_worse,
            'success':bool(superior and both_classes and no_worse)}
