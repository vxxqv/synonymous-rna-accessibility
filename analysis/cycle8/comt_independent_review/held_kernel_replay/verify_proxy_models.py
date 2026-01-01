"""Synthetic verification only. No biological input path is accepted."""
import json, hashlib, random
from pathlib import Path
from fractions import Fraction as F
import proxy_models as m


def main():
    checks=0
    def ok(test):
        nonlocal checks
        if not test: raise AssertionError('Synthetic check '+str(checks+1))
        checks+=1
    x=[-2,-1,0,1,2,3];z=[1,-1,2,-2,0,1];w=[1]*6
    y=[3+2*a+4*b for a,b in zip(x,z)]
    fit=m.fit_linear(x,y,w,z)
    ok((fit['intercept'],fit['rna_slope'],fit['load_slope'])==(3,2,4))
    ok(m.predict_linear(fit,[F(1,2)],[3])==[16])
    null=m.fit_linear(x,[3+2*a for a in x],w,z)
    ok((null['rna_slope'],null['load_slope'])==(2,0))
    constant=m.fit_linear([1]*3,[0,0,9],[1]*3,[2]*3)
    ok(constant['intercept']==3 and m.weighted_median([0,0,9],[1]*3)==0)
    duplicate=m.fit_linear(x,y,w,x)
    ok(duplicate['rna_slope']==0 and duplicate['degeneracy']=='collinear_use_load')
    negative=m.fit_linear(x,[-a for a in x],w)
    ok(negative['rna_slope']==0)
    ok(m.site_weights(['a','a','b'])==[F(1,2),F(1,2),1])
    ok(m.site_mae([0,0,0],[1,1,3],['a','a','b'])==2)
    ok(m.weighted_median([0,2],[1,1])==1)
    p=m.fit_pava([0,0,1,2],[0,2,0,3],[1]*4)
    ok(p['y']==[F(2,3),F(2,3),3])
    ok(m.predict_pava(p,[-1,3])==[F(2,3),3])
    ok(m.range_flags(p,[-1,1,3])==[{'rna_outside':True},{'rna_outside':False},{'rna_outside':True}])
    ok(not m.expansion_criterion(0,[0,1,1,1],125,53)['pass'])
    ok(m.expansion_criterion(9,[10,10,10,10],100,40)['pass'])
    ok(not m.expansion_criterion(9,[10,10,10,10],99,40)['pass'])
    roi=m.fit_roi_medians([0,2,10,12],[1]*4,['ROI_1']*2+['ROI_2']*2)
    ok(m.predict_roi_medians(roi,['ROI_2','ROI_1'])==[11,1])
    try:m.fit_roi_medians([0,1],[1,1],['ROI_1']*2)
    except ValueError:checks+=1
    else:raise AssertionError('Missing ROI accepted')
    # Discarded sequential method is saturated despite identifiable true load.
    sx=[0,1,2,3];sl=[0,1,0,1];sy=[0,2,2,4]
    sp=m.fit_pava(sx,sy,[1]*4)
    residual=[F(a)-b for a,b in zip(sy,m.predict_pava(sp,sx))]
    ok(residual==[0]*4)
    joint=m.fit_linear(sx,sy,[1]*4,sl)
    ok(joint['rna_slope']==1 and joint['load_slope']==1)
    old_prediction=m.predict_pava(sp,[F(1,2)])[0]
    new_prediction=m.predict_linear(joint,[F(1,2)],[0])[0]
    ok(old_prediction==1 and new_prediction==F(1,2))
    # Opposite signs from one shared reference, with no true layer difference.
    d=F(3,2);r=f=h=F(0)
    R,T,L=r-d,d-f,h-r
    ok(R==-T and L==0)
    # Shared RNA noise creates a deterministic apparent load correction.
    r=F(2);d=f=h=F(0)
    R,T,L=r-d,d-f,h-r
    ok(R+L==T==0)
    # True RNA and heavy both equal u: true heavy/RNA load is constant zero.
    # RNA preparation noise e alone yields measured R=u+e and L=-e.
    u=[-1,-1,1,1];noise=[-2,2,-2,2]
    rna=[a+b for a,b in zip(u,noise)];load=[-a for a in noise]
    coupling=m.fit_linear(rna,u,[1]*4,load)
    cp=m.fit_pava(rna,u,[1]*4)
    ok((coupling['rna_slope'],coupling['load_slope'])==(1,1))
    ok(m.site_mae(u,m.predict_linear(coupling,rna,load),list('abcd'))==0)
    ok(m.site_mae(u,m.predict_pava(cp,rna),list('abcd'))==F(1,2))
    rng=random.Random(429)
    for _ in range(100):
        xx=list(range(6)); yy=[rng.randrange(-5,6) for _ in xx]
        ww=[rng.randrange(1,4) for _ in xx]
        got=m.fit_pava(xx,yy,ww)['y']
        # Independent weighted min-max characterization on strictly ordered x.
        ref=[]
        for i in range(6):
            ref.append(max(min(F(sum(ww[k]*yy[k] for k in range(a,b+1)),sum(ww[a:b+1]))
                               for b in range(i,6)) for a in range(i+1)))
        ok(got==ref)
        perm=list(range(6));rng.shuffle(perm)
        pfit=m.fit_pava([xx[i] for i in perm],[yy[i] for i in perm],[ww[i] for i in perm])
        ok(pfit['y']==got)
    for _ in range(80):
        xx=[rng.randrange(-2,3) for _ in range(8)]
        zz=[rng.randrange(-2,3) for _ in range(8)]
        yy=[rng.randrange(-5,6) for _ in range(8)]
        ww=[rng.randrange(1,4) for _ in range(8)]
        fit=m.fit_linear(xx,yy,ww,zz)
        pred=m.predict_linear(fit,xx,zz)
        residual=[F(y)-p for y,p in zip(yy,pred)]
        # KKT conditions independently certify a convex weighted LS optimum.
        ok(sum(t*r for t,r in zip(ww,residual))==0)
        ok(sum(t*z*r for t,z,r in zip(ww,zz,residual))==0)
        cross=sum(t*x*r for t,x,r in zip(ww,xx,residual))
        ok((fit['rna_slope']>0 and cross==0) or (fit['rna_slope']==0 and cross<=0))
        perm=list(range(8));rng.shuffle(perm)
        other=m.fit_linear([xx[i] for i in perm],[yy[i] for i in perm],[ww[i] for i in perm],[zz[i] for i in perm])
        ok([other[k] for k in ('intercept','rna_slope','load_slope')]==[fit[k] for k in ('intercept','rna_slope','load_slope')])
    # Exact rank is for stored floats. Near collinearity must remain visible.
    nf=m.fit_linear([0.,1.,2.,3.],[0,1,2,4],[1]*4,[0.,1.,2.+1e-8,3.])
    ok(nf['centered_rank']==2 and nf['condition_number']>1e8)
    for args in [([],[],[]),([0],[1],[0]),([float('nan')],[1],[1]),([0],[True],[1])]:
        try:m.fit_linear(*args)
        except ValueError:checks+=1
        else:raise AssertionError('Invalid input accepted')
    out={'status':'PASS_SYNTHETIC_ONLY','checks':checks,'source_outcome_rows_read':0,
         'discarded_sequential_pathology':{'true_assessment_target':'1/2','sequential_prediction':'1','joint_prediction':'1/2','training_PAVA_residual':'0'},
         'shared_reference_coupling_controls':['R=-T from gDNA noise','T=R+L from RNA noise without biological load effect','True RNA=heavy=u, RNA noise=e: joint exact MAE0 versus PAVA MAE1/2 with constant true heavy/RNA load'],
         'method_limits':['No trained biological models','Finite numeric input arithmetic is exact, underlying logs are not','No biological uncertainty calibration'],
         'source_sha256':hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest(),
         'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    destination=Path(__file__).with_name('proxy_synthetic_receipt.json')
    if destination.exists():raise FileExistsError('Fresh receipt required')
    destination.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
