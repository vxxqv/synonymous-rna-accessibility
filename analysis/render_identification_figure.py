"""Exact abstract graph examples for the barrier identification results."""
from pathlib import Path
import sys, json, argparse, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'analysis'))
from barrier_identification import barrier,sharp_barrier_bounds
import render_figure_suite as style


def product_graph(g,h):
    return [[w*len(h)+z for w in g[x]]+[x*len(h)+w for w in h[z]] for x in range(len(g)) for z in range(len(h))]


def render(out):
    p3=[[1],[0,2],[1]]
    gh=product_graph(p3,p3)
    g=[1,0,1];h=[0,.5,1]
    vals=[a+b for a in g for b in h]
    route=[0,1,2,5,8,7,6]
    assert barrier(gh,vals,0,6)==0 and min(vals[i] for i in route)==1
    eps=np.linspace(0,.5,41)
    bounds=[sharp_barrier_bounds(p3,[v-e for v in [1,.4,1]],[v+e for v in [1,.4,1]],0,2) for e in eps]
    valley=np.linspace(0,4,41)
    conditional=[];headroom=[]
    g5=[[1],[0,2],[1,3],[2,4],[3]]
    for v in valley:
        hv=[0,-float(v),3]
        conditional.append(barrier(gh,[a+b for a in [2,0,2] for b in hv],0,6))
        headroom.append(barrier(product_graph(g5,p3),[a+b for a in [2,20,0,20,2] for b in hv],0,12))
    cube=[[j for j in range(8) if (i^j).bit_count()==1] for i in range(8)]
    scores=[1,.2,.3,.8,.6,.9,.4,1]
    order=[1,2,4,6,3,5];known={0,7};trace=[]
    for k in range(7):
        lower=[s if i in known else 0 for i,s in enumerate(scores)]
        upper=[s if i in known else 1 for i,s in enumerate(scores)]
        lo,hi=sharp_barrier_bounds(cube,lower,upper,0,7)
        trace.append({'additional_queries':k,'lower':lo,'upper':hi})
        if k<6:known.add(order[k])
    exact=barrier(cube,scores,0,7)
    assert all(r['lower']<=exact<=r['upper'] for r in trace)
    data={'claim':'Sharp score bounds and additive background bypasses on declared finite graphs',
          'scope':'Exact synthetic examples, not biological measurements or evolutionary trajectories',
          'panels':{'a':{'g':g,'h':h,'product_score':vals,'threshold_route':route,'local_barrier':1,'product_barrier':0},
                    'b':{'epsilon':eps.tolist(),'sharp_bounds':bounds,'true_field':[1,.4,1]},
                    'c':{'background_valley':valley.tolist(),'fixed_focal_barrier':conditional,'focal_headroom_barrier':headroom,'formula_without_hypothesis':0},
                    'd':{'cube_scores':scores,'query_order':order,'trace':trace,'true_barrier':exact}}}
    fig=plt.figure(figsize=(174/25.4,169/25.4))
    gs=fig.add_gridspec(2,2,left=.14,right=.96,bottom=.10,top=.90,wspace=.67,hspace=.78)
    ax=fig.add_subplot(gs[0,0],projection='3d')
    ax.view_init(elev=25,azim=-58)
    ax.set_proj_type('ortho')
    for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
        axis.pane.set_alpha(1)
    for u,ns in enumerate(gh):
        for v in ns:
            if u<v:
                ax.plot([u//3,v//3],[h[u%3],h[v%3]],[vals[u],vals[v]],c=style.LIGHT,lw=1)
    ax.plot([0,1,2],[0,0,0],g,color=style.ORANGE,lw=2,marker='o',ms=3)
    ax.plot([u//3 for u in route],[h[u%3] for u in route],[vals[u] for u in route],c=style.TEAL,lw=2.2,marker='o',ms=3)
    ax.set(xlabel='Focal coordinate',ylabel='Background gain',zlabel='Score',xticks=[0,1,2],yticks=[0,1],zticks=[0,1,2])
    ax.tick_params(labelsize=7,pad=0)
    ax.xaxis.labelpad=8;ax.yaxis.labelpad=8;ax.zaxis.labelpad=5
    ax.text2D(0,1.12,'a',transform=ax.transAxes,weight='bold',fontsize=10)
    ax.text2D(.11,1.12,'An additive bypass',transform=ax.transAxes,weight='bold',fontsize=8.5)
    ax=fig.add_subplot(gs[0,1]);style.panel(ax,'b','Sharp simultaneous score bounds')
    b=np.array(bounds)
    ax.fill_between(eps,b[:,0],b[:,1],color='#D7E6E3')
    ax.plot(eps,b[:,0],color=style.TEAL,label='Lower bound')
    ax.plot(eps,b[:,1],color=style.BLUE,label='Upper bound')
    ax.axhline(.6,c=style.GREY,ls=':',label='True barrier')
    ax.set(xlabel=r'Simultaneous score error $\epsilon$',ylabel='Barrier height',ylim=(-.04,1.65),xlim=(0,.5))
    ax.legend(fontsize=7.2,loc='upper left')
    ax=fig.add_subplot(gs[1,0]);style.panel(ax,'c','The background condition matters')
    ax.plot(valley,conditional,color=style.ORANGE,label='Focal path (2, 0, 2)')
    ax.plot(valley,headroom,color=style.BLUE,ls='--',label='Extra focal headroom')
    ax.scatter([0],[0],c=style.TEAL,s=30,zorder=4)
    ax.set(xlabel=r'Background valley depth $v$',ylabel='Actual product barrier',ylim=(-.12,2.2),xlim=(-.1,4.1),yticks=[0,1,2])
    ax.legend(loc='center right',bbox_to_anchor=(1.0,0.25),fontsize=7.0)
    ax=fig.add_subplot(gs[1,1]);style.panel(ax,'d','Missing routes require abstention')
    x=[r['additional_queries'] for r in trace]
    ax.fill_between(x,[r['lower'] for r in trace],[r['upper'] for r in trace],step='post',color='#D7E6E3')
    ax.step(x,[r['lower'] for r in trace],where='post',color=style.TEAL,label='Sharp range')
    ax.step(x,[r['upper'] for r in trace],where='post',color=style.BLUE)
    ax.axhline(exact,c=style.ORANGE,ls=':',label='Complete truth')
    ax.set(xlabel='Additional queried vertices',ylabel='Barrier height',xlim=(0,6),ylim=(-.03,1.35),xticks=[0,2,4,6],yticks=[0,.4,.8,1.2])
    ax.legend(loc='upper right',fontsize=7.2)
    style.SEEN.add('results/programme/identification_v1/checks.json')
    style.SEEN.add('analysis/render_identification_figure.py')
    data['dedicated_renderer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    style.save(fig,'Figure_7',out,data)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    if args.out.resolve().is_relative_to(ROOT):raise ValueError('Figure exports belong outside the computational repository')
    args.out.mkdir(parents=True,exist_ok=True);render(args.out)
