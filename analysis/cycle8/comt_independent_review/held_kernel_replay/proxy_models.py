"""Small pre-outcome model kernel; no source-data reader or outcome execution.

Sufficient statistics are rational arithmetic on supplied finite values. This
does not make logarithms or biological quantities exact. All fitting is local
to the training arguments. No target is accepted by prediction methods.
"""
from fractions import Fraction as F
from collections import Counter
from bisect import bisect_left
import math


def number(x):
    if isinstance(x, bool) or not isinstance(x, (int, float, F)):
        raise ValueError('Expected finite real number')
    if isinstance(x, float) and not math.isfinite(x):
        raise ValueError('Expected finite real number')
    return F(x)


def site_weights(sites):
    if not sites or any(not isinstance(s, str) or not s for s in sites):
        raise ValueError('Nonempty site identifiers required')
    counts = Counter(sites)
    return [F(1, counts[s]) for s in sites]


def validated(x, y, w):
    if not x or len(x) != len(y) or len(x) != len(w):
        raise ValueError('Nonempty aligned training vectors required')
    x, y, w = [list(map(number, z)) for z in (x, y, w)]
    if min(w) <= 0:
        raise ValueError('Positive training weights required')
    return x, y, w


def weighted_median(y, w):
    _, y, w = validated(y, y, w)
    ordered = sorted(zip(y, w))
    half, cum = sum(w)/2, F(0)
    for i, (v, weight) in enumerate(ordered):
        cum += weight
        if cum > half:
            return v
        if cum == half:
            return (v+ordered[i+1][0])/2
    raise AssertionError('Median unreachable')


def fit_pava(x, y, w):
    x, y, w = validated(x, y, w)
    tied = {}
    for a, b, c in zip(x, y, w):
        sw, sy = tied.get(a, (F(0), F(0)))
        tied[a] = sw+c, sy+c*b
    xs = sorted(tied)
    blocks = []
    for i, a in enumerate(xs):
        sw, sy = tied[a]
        blocks.append([i, i, sw, sy])
        while len(blocks) > 1 and blocks[-2][3]/blocks[-2][2] > blocks[-1][3]/blocks[-1][2]:
            right = blocks.pop(); left = blocks.pop()
            blocks.append([left[0], right[1], left[2]+right[2], left[3]+right[3]])
    ys = [None]*len(xs)
    for lo, hi, sw, sy in blocks:
        ys[lo:hi+1] = [sy/sw]*(hi-lo+1)
    return {'kind': 'pava', 'x': xs, 'y': ys, 'blocks': len(blocks),
            'training_range': [min(x), max(x)]}


def fit_roi_medians(y, w, roi):
    _, y, w = validated(y,y,w)
    if len(roi)!=len(y) or any(a not in ('ROI_1','ROI_2') for a in roi):
        raise ValueError('Explicit ROI_1/ROI_2 labels required')
    if set(roi)!={'ROI_1','ROI_2'}:raise ValueError('Each ROI needs training examples')
    return {a:weighted_median([v for v,r in zip(y,roi) if r==a],
                              [v for v,r in zip(w,roi) if r==a])
            for a in ('ROI_1','ROI_2')}


def predict_roi_medians(model, roi):
    if any(a not in model for a in roi):raise ValueError('Unfitted assessment ROI')
    return [model[a] for a in roi]


def predict_pava(model, x):
    out = []
    xs, ys = model['x'], model['y']
    for a in map(number, x):
        i = bisect_left(xs, a)
        if i == 0: out.append(ys[0])
        elif i == len(xs): out.append(ys[-1])
        elif xs[i] == a: out.append(ys[i])
        else: out.append(ys[i-1]+(a-xs[i-1])*(ys[i]-ys[i-1])/(xs[i]-xs[i-1]))
    return out


def fit_linear(x, y, w, load=None):
    x, y, w = validated(x, y, w)
    if load is None: z = [F(0)]*len(x)
    else:
        if len(load) != len(x): raise ValueError('Load alignment')
        z = list(map(number, load))
    sw = sum(w)
    mx, my, mz = [sum(t*v for t, v in zip(w, a))/sw for a in (x, y, z)]
    xc, yc, zc = [[a-m for a in arr] for arr, m in ((x,mx),(y,my),(z,mz))]
    dot = lambda a,b: sum(t*u*v for t,u,v in zip(w,a,b))
    xx, zz, xz, xy, zy = dot(xc,xc), dot(zc,zc), dot(xc,zc), dot(xc,yc), dot(zc,yc)
    det = xx*zz-xz*xz
    if det < 0: raise ArithmeticError('Negative exact Gram determinant')
    if det > 0:
        a, b = (xy*zz-zy*xz)/det, (zy*xx-xy*xz)/det
        rank = 2
        if a < 0: a, b = F(0), zy/zz
        # Stable condition estimate for centered weighted predictor design.
        trace = xx+zz
        ratio = float(det/(trace*trace))
        largest = (1+math.sqrt(max(0.,1-4*ratio)))/2
        cond = None if ratio == 0 else largest/math.sqrt(ratio)
        degeneracy = 'full_rank'
    elif zz > 0:
        a, b, rank, cond, degeneracy = F(0), zy/zz, 1, None, 'collinear_use_load'
    elif xx > 0:
        a, b, rank, cond, degeneracy = max(F(0),xy/xx), F(0), 1, 1., 'load_constant_use_rna'
    else:
        a, b, rank, cond, degeneracy = F(0), F(0), 0, None, 'constant_weighted_mean'
    return {'kind':'linear','intercept':my-a*mx-b*mz,'rna_slope':a,'load_slope':b,
            'centered_rank':rank,'design_rank':rank+1,'condition_number':cond,
            'condition_scope':'centered weighted predictor design, null for singularity or float-conversion underflow',
            'gram_determinant':det,'degeneracy':degeneracy,
            'rna_range':[min(x),max(x)],'load_range':[min(z),max(z)],
            'weighted_response_mean':my}


def predict_linear(model, x, load=None):
    x = list(map(number,x))
    z = [F(0)]*len(x) if load is None else list(map(number,load))
    if len(x) != len(z): raise ValueError('Prediction alignment')
    return [model['intercept']+model['rna_slope']*a+model['load_slope']*b for a,b in zip(x,z)]


def range_flags(model, x, load=None):
    if model['kind']=='pava':
        lo,hi=model['training_range']
        return [{'rna_outside':not lo<=number(a)<=hi} for a in x]
    z = [0]*len(x) if load is None else load
    if len(z)!=len(x):raise ValueError('Range alignment')
    return [{'rna_outside':not model['rna_range'][0]<=number(a)<=model['rna_range'][1],
             'load_outside':not model['load_range'][0]<=number(b)<=model['load_range'][1]} for a,b in zip(x,z)]


def site_mae(observed, predicted, sites):
    if not observed or len(observed)!=len(predicted) or len(observed)!=len(sites):
        raise ValueError('Scoring alignment')
    w = site_weights(sites)
    return sum(t*abs(number(a)-number(b)) for t,a,b in zip(w,observed,predicted))/len(set(sites))


def expansion_criterion(load_mae, baseline_maes, variants, sites):
    if len(baseline_maes)!=4:raise ValueError('Exactly four baseline MAEs required')
    target=number(load_mae); base=list(map(number,baseline_maes))
    if min([target]+base)<0:raise ValueError('MAE must be nonnegative')
    return {'coverage_pass':variants>=100 and sites>=40,
            'relative_comparison_defined':all(b>0 for b in base),
            'gain_pass_each':[b>0 and target<=F(9,10)*b for b in base],
            'pass':variants>=100 and sites>=40 and all(b>0 and target<=F(9,10)*b for b in base)}
