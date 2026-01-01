"""Render frozen Figure 12 panel data without assay access or model fitting."""
import sys
if not __debug__:
    raise RuntimeError('Optimization refused before input access')
sys.dont_write_bytecode = True
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

INK, TEAL, GREY, LIGHT, ORANGE = '#20272B', '#176B68', '#737C84', '#DCE1E4', '#B85E31'
PINS = {'evaluation':'95fb1c15cfe6148c41a543240e1c84db0d54b910cbe974eeec26914202ceb2ac',
        'verification':'9594b34f602488594a855a249a7e34db260de108a7afdf87d1eae3e631840c45',
        'cohort_receipt':'75cef8d4b9954e5f4975e9be077572e6e28217025e56b6d76163b24d12f92318',
        'kinetic_proof':'c9e2a82f62c97949b01f87d541ccdf0ed6ca1e1445f5e9f0cb288471b85cc2b9'}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True)
    parser.add_argument('--expected-data-sha256', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    out = Path(args.out)
    if out.exists():
        raise FileExistsError('Fresh output directory required')
    if sha(args.data) != args.expected_data_sha256:
        raise ValueError('Panel data hash mismatch')
    data = json.loads(Path(args.data).read_text(encoding='utf-8'))
    if data['renderer_sha256'] != sha(__file__) or data['source_sha256'] != PINS:
        raise ValueError('Renderer/source identities differ from panel data')
    if data['primary_success'] is not False or data['analysis_status'] != 'post_exposure_exploratory':
        raise ValueError('Frozen interpretation changed')
    panels = data['panels']
    if len(panels['a']['models']) != 9 or panels['c']['n'] != 5416 or panels['c']['points_excluded'] != 0:
        raise ValueError('Figure coverage changed')
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8.5,
        'axes.labelsize':8.5, 'xtick.labelsize':8.5, 'ytick.labelsize':8.5,
        'text.color':INK, 'axes.labelcolor':INK, 'axes.edgecolor':GREY,
        'xtick.color':INK, 'ytick.color':INK, 'axes.linewidth':.65,
        'xtick.major.width':.65, 'ytick.major.width':.65,
        'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none', 'svg.hashsalt':'Figure12_frozen_v1',
        'savefig.facecolor':'white', 'figure.facecolor':'white'})
    fig = plt.figure(figsize=(174/25.4, 164/25.4))
    fig.text(.025,.975,'Primary model misses the 5% practical criterion', va='top', fontsize=10.5, weight='bold')
    fig.text(.025,.943,'Source-exposed exploratory comparison | 1,112 model-held-out genes', va='top')
    titles = [('a',.025,.89,'All nine frozen models'), ('b',.60,.89,'Matched MAE reductions'),
              ('c',.025,.495,'Coupled contrasts'), ('d',.59,.495,'Synthetic nonidentifiability')]
    for letter,x,y,title in titles:
        fig.text(x,y,letter,fontsize=10,weight='bold',va='center')
        fig.text(x+.034,y,title,fontsize=9.2,weight='bold',va='center')

    # Full nine-model comparison. A narrow numerical scale is explicit.
    axa = fig.add_axes([.31,.61,.15,.245])
    models = panels['a']['models']
    axa.set(xlim=(.33,.49), ylim=(8.7,-.7), xticks=[.34,.40,.48], yticks=[])
    axa.spines[['left','right','top']].set_visible(False)
    axa.set_xlabel('MAE (ln half-life)',labelpad=5)
    axa.grid(axis='x',color=LIGHT,linewidth=.5)
    axa.set_axisbelow(True)
    for i, model in enumerate(models):
        primary = model['role']=='primary'
        sensitivity = model['role']=='sensitivity'
        color = TEAL if primary else GREY
        axa.plot(model['MAE'],i,'o',ms=4.8,mec=color,mew=1,
                 mfc='white' if sensitivity else color,zorder=4)
        fy = fig.transFigure.inverted().transform(axa.transData.transform((.33,i)))[1]
        fig.text(.025,fy,model['label'],va='center',color=color if primary else INK,
                 weight='bold' if primary else 'normal')
        fig.text(.565,fy,f"{model['MAE']:.5f}",ha='right',va='center',color=color,
                 weight='bold' if primary else 'normal')
    sep_y = fig.transFigure.inverted().transform(axa.transData.transform((.33,5.5)))[1]
    fig.add_artist(plt.Line2D([.025,.565],[sep_y,sep_y],transform=fig.transFigure,color=LIGHT,lw=.65))
    fig.text(.025,.525,'Open markers: prespecified sensitivities',color=GREY)

    axb = fig.add_axes([.735,.66,.235,.17])
    axb.set(xlim=(-.15,5.8), ylim=(2.65,-.65), xticks=[0,1,3,5], yticks=[])
    axb.spines[['left','right','top']].set_visible(False)
    axb.set_xlabel('Relative MAE reduction (%)',labelpad=5)
    axb.axvline(5,color=ORANGE,ls=(0,(3,2)),lw=1)
    axb.text(5,-.60,'5%',ha='center',va='bottom',color=ORANGE)
    for i,gain in enumerate(panels['b']['gains']):
        color=TEAL if i==0 else GREY
        axb.plot(gain['percent_gain'],i,'o',ms=5,mec=color,mfc=color if i==0 else 'white',mew=1)
        fy=fig.transFigure.inverted().transform(axb.transData.transform((0,i)))[1]
        fig.text(.60,fy,gain['label'],va='center',color=color)
        axb.annotate(f"{gain['percent_gain']:.3f}%",(gain['percent_gain'],i),
                     xytext=(7,0),textcoords='offset points',va='center',color=color)
    fig.text(.60,.565,'Each gain uses its matched codon baseline.',fontsize=8.5)
    fig.text(.60,.535,'Dashed line: primary practical criterion.',fontsize=8.5)

    axc=fig.add_axes([.11,.125,.345,.315])
    points=panels['c']['all_points']
    xx=[p['lifetime_change_rank_percent'] for p in points]
    yy=[p['rate_contrast_rank_percent'] for p in points]
    axc.scatter(xx,yy,s=2.6,c=TEAL,alpha=.27,edgecolors='none',rasterized=False)
    axc.set(xlim=(-2,102),ylim=(-2,102),xticks=[0,50,100],yticks=[0,50,100],
            xlabel='Lifetime-change rank (%)',ylabel='Reciprocal-contrast rank (%)')
    axc.spines[['top','right']].set_visible(False)
    fig.text(.11,.46,rf"$\rho={panels['c']['rho']:.3f}$; all 5,416 eligible genes",fontsize=8.5)

    axd=fig.add_axes([.59,.12,.385,.325]);axd.set_axis_off()
    axd.set(xlim=(0,1),ylim=(0,1))
    axd.text(0,.97,r'$k_0=b+p,\quad k_1=b+\delta+rp$',va='top',fontsize=10)
    axd.text(0,.79,r'Same rates: $(k_0,k_1)=(3,2)$',va='top')
    axd.text(0,.66,r'Fixed residual activity: $r=1/2$',va='top')
    xs=[.00,.37,.62,.87]
    for x,label in zip(xs,['Setting',r'$b$',r'$p$',r'$\delta$']):
        axd.text(x,.48,label,va='center',weight='bold')
    axd.plot([0,.99],[.405,.405],color=LIGHT,lw=.75,transform=axd.transAxes)
    for y,vals in [(.31,['1','1','2','0']),(.14,['2','0','3',r'$1/2$'])]:
        for j,(x,val) in enumerate(zip(xs,vals)):
            axd.text(x,y,val,va='center',color=TEAL if j==2 else INK,
                     weight='bold' if j==2 else 'normal')
    axd.text(0,-.08,r'Different $p$; identical observed total rates.',va='top')
    fig.text(.025,.022,'Exploratory and source-exposed; no calibrated biological uncertainty.',fontsize=8.5,color=GREY)

    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    texts=[text for text in fig.findobj(matplotlib.text.Text) if text.get_visible() and text.get_text()]
    outer=fig.bbox
    clips=[]
    for text in texts:
        bounds=text.get_window_extent(renderer)
        if bounds.x0<outer.x0-.5 or bounds.y0<outer.y0-.5 or bounds.x1>outer.x1+.5 or bounds.y1>outer.y1+.5:
            clips.append(text.get_text())
    if clips:
        raise ValueError('Text outside page: '+repr(clips))
    minimum=min(t.get_fontsize() for t in texts)
    if minimum<8.5:
        raise ValueError('Text below 8.5 point minimum')
    out.mkdir(parents=True)
    metadata={'Creator':'','Producer':'','CreationDate':None,'ModDate':None}
    fig.savefig(out/'Figure_12.pdf',metadata=metadata)
    fig.savefig(out/'Figure_12.svg',metadata={'Creator':'','Date':None})
    fig.savefig(out/'Figure_12.png',dpi=610,metadata={'Software':''})
    fig.savefig(out/'Figure_12_preview.png',dpi=180,metadata={'Software':''})
    receipt={'schema':'Figure12_render_receipt_v1','data_sha256':sha(args.data),'renderer_sha256':sha(__file__),
       'physical_dimensions_mm':[174,164],'minimum_text_pt':minimum,'visible_text_count':len(texts),
       'text_outside_page':clips,'models_plotted':len(models),'kinetic_points_plotted':len(xx),
       'points_excluded':0,'x_plot_range_rank_percent':[-2,102],'y_plot_range_rank_percent':[-2,102],
       'manual_visual_review':'Pending separate review report',
       'asset_sha256':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
    (out/'render_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    plt.close(fig)
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    main()
