"""Publication Figure 8: exact synthetic correlated uncertainty examples.

All numerical marks derive from the exact affine API, compared with its held
release fixtures. No statistical calibration, sampling or biological data.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import itertools
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis' / 'programme'))
from joint_uncertainty_affine import analyze, graph_from_edges, score_at, widest_path

INK, TEAL, GREY, LIGHT = '#20272B', '#176B68', '#737C84', '#DCE1E4'
INPUTS = ['analysis/programme/joint_uncertainty_affine.py',
          'results/programme/joint_uncertainty_release/joint_uncertainty_results.json',
          'results/programme/joint_uncertainty_release/manifest.json']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_json(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {k: exact_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [exact_json(v) for v in value]
    return value


def panel(ax, letter, title):
    ax.text(-.015, 1.15, letter, transform=ax.transAxes, fontsize=10,
            weight='bold', va='bottom')
    ax.text(.065, 1.15, title, transform=ax.transAxes, fontsize=8.5,
            weight='bold', va='bottom')


def render(out, qa):
    if not __debug__:
        raise RuntimeError('Figure integrity checks require normal Python execution')
    before = {p: sha(ROOT / p) for p in INPUTS}
    coefficients = [(1, 0), (0, 1), (0, -1), (1, 0)]
    shifted = [(4, 2), (3, 3), (3, 1), (4, 2)]
    edges = [(0, 1), (1, 3), (0, 2), (2, 3)]
    base = analyze(4, edges, coefficients, 0, 3, -1, 1)
    offset = analyze(4, edges, shifted, 0, 3, -1, 1)
    release = json.loads((ROOT / INPUTS[1]).read_text())
    assert exact_json(base) == release['diamond']
    assert exact_json(offset) == release['shared_affine_offset']
    theta = [Fraction(k, 100) for k in range(-100, 101)]
    graph = graph_from_edges(4, edges)
    curves = []
    for t in theta:
        f = score_at(coefficients, t)
        f2 = score_at(shifted, t)
        capacity, path = widest_path(graph, f, 0, 3)
        capacity2, _ = widest_path(graph, f2, 0, 3)
        barrier = min(f[0], f[3]) - capacity
        barrier2 = min(f2[0], f2[3]) - capacity2
        assert capacity == abs(t) and barrier == 1-abs(t)
        assert capacity2 == capacity+3+2*t and barrier2 == barrier
        curves.append({'theta': t, 'capacity': capacity, 'barrier': barrier,
                       'shifted_capacity': capacity2, 'shifted_barrier': barrier2,
                       'optimal_path': path})
    x = np.array([float(t) for t in theta])
    capacity = np.array([float(v['capacity']) for v in curves])
    barrier = np.array([float(v['barrier']) for v in curves])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .65, 'text.color': INK, 'axes.labelcolor': INK,
        'xtick.color': INK, 'ytick.color': INK, 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    fig, axs = plt.subplots(3, 2, figsize=(174/25.4, 153/25.4))
    fig.subplots_adjust(left=.13, right=.97, bottom=.085, top=.91,
                        wspace=.62, hspace=1.05)

    ax = axs[0, 0]; panel(ax, 'a', 'Two competing routes')
    pos = {0: (.08, .5), 1: (.5, .87), 2: (.5, .13), 3: (.92, .5)}
    for u, v in edges:
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=GREY, lw=1.05, zorder=1)
    for v, (xx, yy) in pos.items():
        ax.scatter(xx, yy, s=23, facecolor='white', edgecolor=INK, lw=.9, zorder=3)
    for xx, yy, text, ha in [(-.02,.5,r'$s:1$','right'),
                            (1.02,.5,r'$t:1$','left'),
                            (.5,1.04,r'$u:\theta$','center'),
                            (.5,-.04,r'$v:-\theta$','center')]:
        ax.text(xx, yy, text, ha=ha, va='center', fontsize=8.5)
    ax.text(.5,-.35,r'$-1\leq\theta\leq1$',ha='center',va='center')
    ax.set(xlim=(-.2,1.2), ylim=(-.2,1.2)); ax.axis('off')

    ax = axs[0, 1]; panel(ax, 'b', 'Joint line and marginal box')
    ax.add_patch(Rectangle((-1,-1),2,2,facecolor='#F2F4F4',
                           edgecolor=GREY,lw=.9,ls='--'))
    ax.plot(x, -x, c=TEAL, lw=1.8)
    ax.scatter([-1,1], [1,-1], c=TEAL, s=16, zorder=4)
    ax.scatter([-1], [-1], marker='s', facecolor='white', edgecolor=GREY, s=22)
    ax.text(.70,.75,'Box',transform=ax.transAxes,color=GREY)
    ax.text(.12,.22,'Joint line',transform=ax.transAxes,color=TEAL)
    ax.set(xlim=(-1.2,1.2),ylim=(-1.2,1.2),xticks=[-1,0,1],yticks=[-1,0,1],
           xlabel='Score at u',ylabel='Score at v')

    ax = axs[1, 0]; panel(ax, 'c', 'The largest barrier is interior')
    ax.plot(x, capacity, c=TEAL, lw=1.7)
    ax.plot(x, barrier, c=INK, lw=1.3, ls='--')
    ax.scatter([-1,0,1], [1,0,1], c=TEAL, s=14, zorder=3)
    ax.scatter([-1,0,1], [0,1,0], c=INK, marker='s', s=12, zorder=3)
    ax.legend([Line2D([],[],color=TEAL,lw=1.7),
               Line2D([],[],color=INK,lw=1.3,ls='--')],
              [r'$C=|\theta|$',r'$B=1-|\theta|$'],
              loc='lower left',bbox_to_anchor=(0,.91),ncol=2,frameon=False,
              fontsize=8,handlelength=1.2,columnspacing=1.2,borderaxespad=0)
    ax.set(xlim=(-1.1,1.1),ylim=(-.08,1.18),xticks=[-1,0,1],yticks=[0,.5,1],
           xlabel=r'Joint parameter $\theta$',ylabel='Capacity / barrier')

    ax = axs[1, 1]; panel(ax, 'd', 'Changing the order of decisions')
    yy = np.array([2,1,0])
    joint_values = [base[k] for k in ('fixed','guaranteed','possible')]
    box_values = [base[k] for k in ('fixed','box_guaranteed','box_possible')]
    ax.scatter([float(v) for v in joint_values], yy+.1, c=TEAL, s=21, zorder=3)
    ax.scatter([float(v) for v in box_values], yy-.1, facecolor='white',
               edgecolor=GREY, marker='s', s=21, zorder=3)
    for value, y in zip(joint_values, yy+.1):
        ax.text(float(value)+.09,y,str(value).replace('-','−'),va='center',color=TEAL)
    ax.set(yticks=yy,yticklabels=['Fixed','Guaranteed','Possible'],
           xticks=[-1,0,1],xlim=(-1.25,1.42),ylim=(-.4,2.8),
           xlabel='Absolute route capacity')
    ax.legend([Line2D([],[],marker='o',ls='',color=TEAL,ms=4),
               Line2D([],[],marker='s',ls='',mfc='white',mec=GREY,ms=4)],
              ['Joint','Box'],ncol=2,loc='upper left',fontsize=8,
              frameon=False,handlelength=.6,columnspacing=.8,borderaxespad=0)

    ax = axs[2, 0]; panel(ax, 'e', 'Box bounds and endpoint-only failure')
    ranges = [base['barrier_range'],base['box_barrier_range'],(0,0)]
    for y, (lo, hi), color, ls in zip([2,1,0],ranges,[TEAL,GREY,INK],['-','--',':']):
        ax.plot([float(lo),float(hi)],[y,y],c=color,lw=2,ls=ls)
        ax.plot([float(lo),float(hi)],[y,y],'|',c=color,ms=7,mew=1.2)
        ax.text(float(hi)+.12,y,('[0, '+str(hi)+']') if hi else '{0}',va='center',color=color)
    ax.set(yticks=[2,1,0],yticklabels=['Joint set','Full box','θ = ±1 only'],
           xticks=[0,1,2],xlim=(-.1,2.75),ylim=(-.4,2.6),
           xlabel='Endpoint-relative barrier')

    ax = axs[2, 1]; panel(ax, 'f', 'Shared shifts preserve B')
    ax.axis('off'); ax.set(xlim=(0,1),ylim=(0,1))
    ax.text(.52,.94,'Original',ha='center')
    ax.text(.87,.94,r'$+\,3+2\theta$',ha='center')
    ax.plot([0,1],[.84,.84],c=LIGHT,lw=.8)
    rows = [('Fixed','fixed'),('Guaranteed','guaranteed'),
            ('Possible','possible'),('Max B','barrier_max')]
    for y,(label,key) in zip([.71,.51,.31,.11],rows):
        ax.text(0,y,label,va='center')
        for xx,result in [(.52,base),(.87,offset)]:
            value = result[key] if key!='barrier_max' else result['barrier_range'][1]
            ax.text(xx,y,str(value).replace('-','−'),ha='center',va='center')

    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(matplotlib.text.Text)
           if t.get_visible() and t.get_text().strip()]
    clipping=[]; overlaps=[]
    for t in texts:
        b=t.get_window_extent(renderer)
        if not fig.bbox.contains(b.x0,b.y0) or not fig.bbox.contains(b.x1,b.y1):
            clipping.append(t.get_text())
        assert t.get_fontsize() >= 8
    for a,b in itertools.combinations(texts,2):
        ba,bb=a.get_window_extent(renderer),b.get_window_extent(renderer)
        if min(ba.x1,bb.x1)-max(ba.x0,bb.x0)>.5 and min(ba.y1,bb.y1)-max(ba.y0,bb.y0)>.5:
            overlaps.append([a.get_text(),b.get_text()])
    if clipping or overlaps:
        raise RuntimeError({'clipping':clipping,'text_overlaps':overlaps})
    out.mkdir(parents=True,exist_ok=True);qa.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','svg','eps','png']:
        fig.savefig(out/f'Figure_8.{ext}',dpi=610 if ext=='png' else 600)
    fig.savefig(qa/'final_size.png',dpi=180)
    plt.close(fig)
    from PIL import Image
    Image.open(qa/'final_size.png').convert('L').save(qa/'grayscale.png')
    after={p:sha(ROOT/p) for p in INPUTS}
    assert before==after
    data={'claim':'Correlation and decision order change absolute route guarantees; endpoint checks can miss the maximum relative barrier',
          'scope':'Exact synthetic graph, deterministic closed uncertainty sets, no confidence intervals',
          'graph':{'vertices':['s','u','v','t'],'edges':edges,'coefficients':coefficients,'parameter_interval':[-1,1]},
          'diamond':exact_json(base),'shared_affine_offset':exact_json(offset),
          'curves':exact_json(curves),
          'panels':{'a':'Score-labeled diamond','b':{'joint_constraint':'f_v = -f_u','intermediate_marginals':[-1,1]},
                    'c':{'capacity':'abs(theta)','barrier':'1-abs(theta)'},
                    'd':{'order':['fixed','guaranteed','possible'],'joint':exact_json(joint_values),'box':exact_json(box_values)},
                    'e':{'joint_barrier_range':exact_json(base['barrier_range']),'box_barrier_range':exact_json(base['box_barrier_range']),'interval_endpoints_only':[0,0]},
                    'f':{'common_offset':'3+2theta','barrier_invariant':True}},
          'export':{'width_mm':174,'height_mm':153,'minimum_text_pt':8,'raster_dpi':610,'clipped_text':clipping,'text_overlaps':overlaps},
          'provenance':{'input_hashes_before':before,'input_hashes_after':after,'renderer_sha256':sha(Path(__file__)),
                        'matplotlib':matplotlib.__version__,'numpy':np.__version__,
                        'release_fixture_equal':True,'exact_curve_checks':len(curves)}}
    destination=ROOT/'results'/'figure_data'/'Figure_8.json'
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (qa/'geometry_qa.json').write_text(json.dumps(data['export'],indent=2)+'\n')
    print(json.dumps({'figure':str(out/'Figure_8.pdf'),'data':str(destination),'exact_points':len(curves),'clipping':clipping,'overlaps':overlaps}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--qa',type=Path,required=True);args=p.parse_args()
    if args.out.resolve().is_relative_to(ROOT):
        raise ValueError('Figure exports must be outside the computational repository')
    render(args.out,args.qa)
