"""Render the frozen multiclass measurement comparison from recorded results."""
from pathlib import Path
from itertools import combinations
import argparse
import hashlib
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RESULT='results/programme/gill_multiclass_v1/validate/validation.json'
PROTO='analysis/programme/retrospective/gill_multiclass_v1/protocol.json'
EXPECTED={RESULT:'cdbdd83d6861c54d1cd5094990cbeb306803349c63bdd0519caf391676c819f6',
          PROTO:'4ab4e85c14203c0ab8c983e8433c766f1566f373a335fd3fda4fb3981371c1aa'}
INK='#20272B';TEAL='#176B68';PURPLE='#735C91';ORANGE='#B86D2D';GREY='#6B7580';LIGHT='#DCE1E4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--qa',type=Path,required=True)
    args=ap.parse_args()
    out=args.outdir.resolve();qa=args.qa.resolve()
    if out.is_relative_to(ROOT):raise ValueError('Figure exports must stay outside the computational repository')
    out.mkdir(parents=True,exist_ok=True);qa.mkdir(parents=True,exist_ok=True)
    for name,digest in EXPECTED.items():assert sha(ROOT/name)==digest
    p=json.loads((ROOT/PROTO).read_text());d=json.loads((ROOT/RESULT).read_text())
    primary=d['primary_result'];loss=primary['losses']
    assert primary['class_count']==1052 and primary['status']=='target_failed'
    assert abs(primary['cosine']-.07787615736918294)<1e-14
    terms=primary['primary_decomposition_aggregate_terms']
    assert abs(math.fsum(terms[:4])-terms[4])<1e-10
    influences=primary['leave_one_stratum_out']
    assert len(influences)==154 and all(x['remaining_cosine'] is not None for x in influences)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
       'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,
       'axes.edgecolor':GREY,'axes.linewidth':.6,'text.color':INK,'axes.labelcolor':INK,
       'xtick.color':INK,'ytick.color':INK,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
    fig=plt.figure(figsize=(174/25.4,170/25.4),facecolor='white')
    positions=[(.15,.755,.31,.145),(.63,.755,.31,.145),(.15,.455,.31,.17),
               (.63,.455,.31,.17),(.15,.15,.31,.17),(.63,.15,.31,.17)]
    axes=[fig.add_axes(x) for x in positions]
    titles=['Frozen protein-preserving groups','Residual repeatability fails',
            'Count agreement is insufficient','Large terms nearly cancel',
            'Fixed predictions lose to baselines','Failure is not one stratum']
    for ax,letter,title in zip(axes,'abcdef',titles):
        ax.text(-.12,1.23,letter,transform=ax.transAxes,fontsize=10,weight='bold')
        ax.text(0,1.23,title,transform=ax.transAxes,fontsize=8.5,weight='bold')
        ax.tick_params(length=3,width=.6)
    a,b,c,equ,e,f=axes
    # Eligibility is frozen from sequence metadata. Bars are group counts.
    a.barh([1,0],[p['development_class_count'],p['primary_class_count']],height=.46,color=[TEAL,PURPLE])
    a.set_yticks([1,0],['Development','Validation'])
    a.set_xlim(0,1450);a.set_xticks([0,500,1000]);a.set_xlabel('Eligible groups')
    a.set_ylim(-.6,1.7)
    for y,n,vertices in [(1,1124,20028),(0,1052,18310)]:
        a.text(n+30,y,f'{n:,}',va='center',fontsize=8)
        a.text(40,y+.34,f'{vertices:,} genotypes',fontsize=8,color=GREY)
    a.text(0,1.08,'38,338 genotypes; one reporter library',transform=a.transAxes,fontsize=8,color=GREY)
    # Prespecified pseudocount settings, all retained.
    for color,partition,offset in [(TEAL,'development',.10),(PURPLE,'validation',-.10)]:
        values=[d['all_summaries'][pc][partition]['cosine'] for pc in ['0.25','0.5','1.0']]
        b.plot(values,np.arange(3)+offset,'o',color=color,markersize=4)
    b.axvline(.5,color=ORANGE,lw=1,ls=(0,(3,2)))
    b.text(.50,2.49,'Target 0.5',ha='center',color=ORANGE,fontsize=8)
    b.set_xlim(0,.6);b.set_ylim(-.45,2.85);b.set_xticks([0,.25,.5])
    b.set_yticks([0,1,2],['0.25','0.5','1']);b.set_ylabel('Pseudocount')
    b.set_xlabel('Pooled residual cosine')
    # Three RNA and three DNA pairwise comparisons, no confidence intervals.
    for y,values,color in [(2,primary['RNA_residual_pooled_cosines'],TEAL),
                            (1,primary['DNA_residual_pooled_cosines'],GREY)]:
        for value,off in zip(values,[-.12,0,.12]):c.plot(value,y+off,'o',color=color,ms=4)
    c.plot(primary['cosine'],0,'o',color=PURPLE,ms=5)
    c.set_yticks([0,1,2],['Normalized','DNA only','RNA only'])
    c.set_xlim(0,1);c.set_xticks([0,.5,1]);c.set_ylim(-.5,2.5)
    c.set_xlabel('Pooled codon-residual cosine')
    # Deterministic signed inner-product identity, not causal attribution.
    x=np.arange(5);equ.bar(x,terms,color=[TEAL,TEAL,ORANGE,ORANGE,PURPLE],width=.65)
    equ.axhline(0,color=GREY,lw=.6)
    for i,value in enumerate(terms):
        equ.text(i,value+(16 if value>=0 else -16),f'{value:.1f}',ha='center',
                 va='bottom' if value>=0 else 'top',fontsize=8)
    equ.set_xticks(x,['RNA\nRNA','DNA\nDNA','RNA\nDNA','DNA\nRNA','Total'])
    equ.set_ylim(-450,450);equ.set_yticks([-400,0,400]);equ.set_ylabel('Weighted cross product')
    # Shape comparisons are centered within class, including validation centering.
    keys=['zero_residual_mse','residual_prediction_mse','constant_prediction_mse',
          'codon_prediction_mse','whole_prediction_mse']
    rmses=[math.sqrt(loss[k]) for k in keys];xx=[0,1.2,3,4.2,5.4]
    e.bar(xx,rmses,width=.70,color=[GREY,PURPLE,GREY,TEAL,PURPLE])
    for x,y in zip(xx,rmses):e.text(x,y+.012,f'{y:.3f}',ha='center',fontsize=8)
    e.set_xticks(xx,['Zero','Fixed','Mean','Add.','Full'])
    e.set_ylim(0,.52);e.set_yticks([0,.25,.5]);e.set_ylabel('Equal-group RMSE')
    e.text(.15,-.30,'Residual',ha='center',transform=e.transAxes,fontsize=8,color=GREY)
    e.text(.73,-.30,'Centered field',ha='center',transform=e.transAxes,fontsize=8,color=GREY)
    # Complete fixed influence diagnostic. No intervals are inferential here.
    influence_values=[x['remaining_cosine'] for x in influences]
    f.scatter(np.arange(1,len(influences)+1),influence_values,s=7,color='#9685AB',linewidths=0)
    f.axhline(primary['cosine'],lw=.8,color=INK,ls=(0,(3,2)))
    f.set_xlim(0,155);f.set_xticks([1,75,154]);f.set_ylim(.05,.11);f.set_yticks([.05,.08,.11])
    f.set_xlabel('Omitted amino-acid stratum');f.set_ylabel('Residual cosine')
    f.text(.5,1.04,'All 154 omissions: 0.068 to 0.092',ha='center',transform=f.transAxes,fontsize=8,color=GREY)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=TEAL,label='Development'),
                        Line2D([],[],marker='o',ls='',color=PURPLE,label='Validation')],
               loc='upper center',bbox_to_anchor=(.53,.995),ncol=2,frameon=False,fontsize=8)
    fig.text(.06,.026,'Fixed cohort and settings. Descriptive library agreement; no calibrated biological uncertainty.',fontsize=8)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    clipping=[];overlaps=[]
    for t in texts:
        box=t.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1>fig.bbox.x1+1 or box.y1>fig.bbox.y1+1:clipping.append(t.get_text())
        if t.get_fontsize()<8:raise ValueError('Ordinary font below8pt')
    for t,u in combinations(texts,2):
        x=t.get_window_extent(renderer);y=u.get_window_extent(renderer)
        if min(x.x1,y.x1)-max(x.x0,y.x0)>.5 and min(x.y1,y.y1)-max(x.y0,y.y0)>.5:
            overlaps.append([t.get_text(),u.get_text()])
    audit={'clipping':clipping,'overlaps':overlaps}
    (qa/'layout.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    fig.savefig(qa/'review.png',dpi=180)
    if clipping or overlaps:raise ValueError(json.dumps(audit))
    for ext in ['pdf','svg','eps','png']:
        metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else {'Date':None} if ext=='svg' else None
        fig.savefig(out/('Figure_10.'+ext),dpi=610,metadata=metadata)
    for name,digest in EXPECTED.items():assert sha(ROOT/name)==digest
    panel={'source_hashes':EXPECTED,'renderer_sha256':sha(__file__),'primary':primary,
        'all_cosines':{pc:{part:s['cosine'] for part,s in v.items()} for pc,v in d['all_summaries'].items()},
        'plotted_losses':dict(zip(keys,rmses)),'cohort':{'development_groups':1124,'validation_groups':1052,
        'development_vertices':20028,'validation_vertices':18310},
        'export':{'width_mm':174,'height_mm':170,'minimum_ordinary_font_pt':8,**audit,
                  'artifacts':{ext:sha(out/('Figure_10.'+ext)) for ext in ['pdf','svg','eps','png']}},
        'scope':'Frozen processed-library results. No confidence intervals, biological replicate independence or causal normalization fractions.'}
    (ROOT/'results/figure_data/Figure_10.json').write_text(json.dumps(panel,indent=2)+'\n',encoding='utf-8')
    (qa/'receipt.json').write_text(json.dumps({'renderer_sha256':sha(__file__),'source_hashes':EXPECTED,'export':panel['export']},indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'figure':str(out/'Figure_10.pdf'),'clipping':clipping,'overlaps':overlaps}))

if __name__=='__main__':main()
