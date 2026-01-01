"""Synthetic-only adversarial checks of existing section 3. No outcome reader."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from pathlib import Path
import hashlib
import json
import sys

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def implied_interval(c, d, interaction):
    return c[0]-d[1]-interaction, c[1]-d[0]+interaction

def construct(alpha, c, d, interaction):
    """Select a feasible difference, then explicit C,D and interaction values."""
    lo, hi = c[0]-d[1], c[1]-d[0]
    difference = min(hi, max(lo, alpha))
    marker = max(d[0], c[0]-difference)
    combined = marker+difference
    eta = difference-alpha
    require(c[0] <= combined <= c[1], 'Combined value outside interval')
    require(d[0] <= marker <= d[1], 'Marker value outside interval')
    require(abs(eta) <= interaction, 'Interaction cap violated')
    require(combined-marker == alpha+eta, 'Contrast reconstruction failed')
    return {'combined':combined, 'marker_only':marker, 'interaction':eta}

def main():
    if not __debug__:
        raise RuntimeError('Optimized Python refused')
    path = Path(__file__).with_name('receipt.json')
    require(not path.exists(), 'Refusing overwrite')
    intervals = list(combinations_with_replacement([F(-1),F(0),F(1)], 2))
    alphas = [F(i,4) for i in range(-16,17)]
    grid = [F(i,4) for i in range(-4,5)]
    options=[]
    single_checks = constructions = 0
    for c,d,cap in product(intervals,intervals,[F(0),F(1,2)]):
        lower,upper=implied_interval(c,d,cap)
        # Independent finite oracle: all C,D grid values, then permitted alpha.
        points={(cc,dd) for cc,dd in product(grid,repeat=2) if c[0]<=cc<=c[1] and d[0]<=dd<=d[1]}
        oracle={a for a in alphas if any(abs(cc-dd-a)<=cap for cc,dd in points)}
        analytic={a for a in alphas if lower<=a<=upper}
        require(oracle==analytic,'Single-marker interval formula mismatch')
        for a in analytic:
            construct(a,c,d,cap)
            constructions+=1
        single_checks+=len(alphas)
        options.append((c,d,cap,oracle,(lower,upper)))
    pair_checks=0
    for a,b in product(options,repeat=2):
        lo=max(a[4][0],b[4][0]);hi=min(a[4][1],b[4][1])
        require(a[3]&b[3]=={v for v in alphas if lo<=v<=hi},'Shared-alpha intersection mismatch')
        if lo<=hi:
            for alpha in (lo,hi):
                construct(alpha,*a[:3]);construct(alpha,*b[:3]);constructions+=2
        pair_checks+=1
    # Arbitrarily strong marker-only calibration and identical observed states
    # still leave opposite isolated effects when interactions are unrestricted.
    worlds=[]
    for alpha in [F(-10),F(10)]:
        mu=F(0);effects={'e':alpha,'a':F(0),'b':F(0)}
        observed={('a','e'):F(2),('b','e'):F(3),('a','b'):F(4)}
        eta={pair:value-mu-effects[pair[0]]-effects[pair[1]] for pair,value in observed.items()}
        predicted={pair:mu+effects[pair[0]]+effects[pair[1]]+eta[pair] for pair in observed}
        require(predicted==observed,'Unordered physical-pair alias consistency')
        require(mu+effects['a']==0 and mu+effects['b']==0,'Marker-only controls changed')
        worlds.append({'alpha_log2':alpha,'WT_ratio':1,'marker_only_ratios':[1,1],
                       'edit_only_ratio':F(2)**alpha.numerator,
                       'combined_log2':{'/'.join(p):v for p,v in predicted.items()},
                       'interactions_log2':{'/'.join(p):v for p,v in eta.items()}})
    require(worlds[0]['combined_log2']==worlds[1]['combined_log2'],'Observational equivalence failed')
    # Known positive-raw WT: main-marker cap alone must respect all raw cells.
    mu=F(1);combined=F(1);alpha=F(2);B=F(2);I=F(0)
    nominal=(combined-mu-B-I,combined-mu+B+I)
    beta=combined-mu-alpha
    require(nominal[0]<=alpha<=nominal[1] and mu+alpha>0 and abs(beta)<=B and mu+beta<0,
            'Raw-positivity counterexample not demonstrated')
    # An interval intersection may identify a nonempty bounded set under an
    # actual absolute interaction cap; both endpoints are attainable.
    positive=[]
    for alpha in [F(1),F(2)]:
        witnesses=[construct(alpha,(F(2),F(2)),(F(1),F(1)),F(1)),
                   construct(alpha,(F(3),F(3)),(F(1),F(1)),F(1))]
        positive.append({'alpha':alpha,'witnesses':witnesses})
    # Latent no-interaction plus monotone nonlinearity need not be additive on
    # the observed log scale: h(z)=2^(z*z) for the used nonnegative levels.
    latent={'00':0,'e0':1,'0m':1,'em':2}
    raw={k:2**(z*z) for k,z in latent.items()}
    observed_log2={k:z*z for k,z in latent.items()}
    require(latent['em']-latent['e0']-latent['0m']+latent['00']==0,'Latent additive fixture')
    require(observed_log2['em']-observed_log2['e0']-observed_log2['0m']+observed_log2['00']==2,
            'Unknown monotone readout counterexample')
    result={'status':'passed','scale':'Exact rational synthetic examples; no empirical bounds',
            'single_marker_alpha_membership_checks':single_checks,'two_marker_systems':pair_checks,
            'constructive_witnesses':constructions,'opposite_effect_observationally_identical_worlds':worlds,
            'raw_main_bound_domain_counterexample':{'mu':mu,'combined':combined,'B':B,'I':I,'nominal_alpha_interval':nominal,
                'admitted_alpha':F(2),'implied_beta':beta,'invalid_marker_only':mu+beta,
                'correct_sharp_alpha_set':'(-1,1) under these strictly positive raw-cell assumptions'},
            'bounded_interaction_positive_example':{'sharp_alpha_set':['1','2'],'endpoint_witnesses':positive},
            'monotone_latent_vs_observed_log':{'latent':latent,'raw':raw,'observed_log2':observed_log2,'observed_log_interaction':2},
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Finite corroboration of independently reconstructed algebra, not a universal proof or real-data calibration.'}
    path.write_text(json.dumps(result,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','single_marker_checks':single_checks,'two_marker_systems':pair_checks,'constructive_witnesses':constructions}))

if __name__=='__main__':
    main()
