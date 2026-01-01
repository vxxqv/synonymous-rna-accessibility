"""Render only frozen FAS outcomes; no model fitting or outcome selection."""
from pathlib import Path
from fractions import Fraction as F
import argparse, gzip, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

parser=argparse.ArgumentParser()
parser.add_argument('--result',required=True,type=Path)
parser.add_argument('--out',required=True,type=Path)
args=parser.parse_args()
raw=args.result.read_bytes()
if args.result.suffix=='.gz':raw=gzip.decompress(raw)
d=json.loads(raw)
args.out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,
    'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':0.7})
fig=plt.figure(figsize=(8.5,9.2),facecolor='white')
gs=fig.add_gridspec(2,2,left=.085,right=.97,bottom=.20,top=.78,
    height_ratios=[1,1.13],width_ratios=[1.16,1],hspace=.55,wspace=.53)
teal='#007F86'; navy='#18354B'; red='#BD4C37'; grey='#95A2AD'
fig.text(.085,.965,'FAS: measured ordering and conditional accessibility',fontsize=14,weight='bold',color=navy)
fig.text(.085,.932,'Native-protein class within the published design | included-RNA enrichment',fontsize=9.5,color=navy)
for x,value,label in [(.085,'72 / 72','finite genotype envelopes'),(.405,'74','forced score orders'),(.71,'0 / 2,556','envelope-certified barrier pairs')]:
    fig.text(x,.88,value,fontsize=18,weight='bold',color=teal)
    fig.text(x,.853,label,fontsize=9,color=navy)

ax=fig.add_subplot(gs[0,:])
rows=d['vertices_data'];x=np.arange(72)
lo=np.array([float(F(v['lower'])) for v in rows]);hi=np.array([float(F(v['upper'])) for v in rows]);pt=np.array([float(F(v['point'])) for v in rows])
ax.vlines(x,lo,hi,color=grey,lw=1.35,zorder=1,label='Range across 27 dependent scenarios')
ax.scatter(x,pt,color=teal,s=12,zorder=3,label='Corrected source-formula point')
ax.axhline(1,color=navy,lw=.65,ls=':',zorder=0)
ax.scatter([46],[1],facecolors='white',edgecolors=navy,s=42,zorder=4)
ax.annotate('Native reference',xy=(46,1),xytext=(50,1.3),fontsize=9,color=navy,arrowprops={'arrowstyle':'-','lw':.6,'color':navy})
ax.set(xlim=(-1.5,72),xlabel='Fixed genotype index (all 72 variants)',ylabel='Relative included-RNA\nenrichment',ylim=(0,max(hi)*1.12))
ax.set_xticks([0,12,24,36,48,60,71]);ax.yaxis.set_major_locator(MaxNLocator(5))
ax.grid(axis='y',color='#E5E9EC',lw=.6)
ax.text(-.075,1.10,'a',transform=ax.transAxes,fontsize=14,weight='bold',color=navy)
ax.set_title('Every declared genotype is retained',loc='left',pad=10,color=navy)
ax.legend(loc='upper left',frameon=False,fontsize=9,handlelength=1.5)

ax=fig.add_subplot(gs[1,0]);pairs=d['primary_codon_null']['dual']
for k,w in enumerate(pairs):
    i,j=w['lower_vertex'],w['higher_vertex'];a=float(F(rows[i]['upper']));b=float(F(rows[j]['lower']))
    y=4-k
    ax.plot([a,b],[y,y],color=red,lw=2.4,zorder=1)
    ax.scatter([a],[y],s=35,color=navy,marker='s',zorder=3)
    ax.scatter([b],[y],s=35,color=teal,marker='o',zorder=3)
    ax.annotate(f'{100*(b/a-1):.2f}%',(max(a,b),y),xytext=(7,0),textcoords='offset points',va='center',fontsize=9,color=red)
ax.set_yticks(range(5),[f"{w['lower_vertex']} < {w['higher_vertex']}" for w in reversed(pairs)])
ax.set(xlim=(.50,1.27),ylim=(-.6,4.6),xlabel='Upper score of lower variant;\nlower score of higher variant')
ax.grid(axis='x',color='#E5E9EC',lw=.6)
ax.set_title('Codon-additive scores with a common\nmonotone response are excluded',loc='left',pad=14,fontsize=10,color=navy)
ax.text(-.20,1.10,'b',transform=ax.transAxes,fontsize=14,weight='bold',color=navy)
fig.text(.085,.065,'Each feature-row difference has weight 1/5.\nTheir sum is zero; all required margins are positive.\nPercentages show separation ratios, not confidence.',fontsize=9,linespacing=1.4,color=navy)

ax=fig.add_subplot(gs[1,1]);arr=np.array([z['bounded_pair_count'] for z in d['scenario_results'][:27]]).reshape(9,3)
im=ax.imshow(arr,cmap='Blues',vmin=0,vmax=320,aspect='auto')
for r in range(9):
    for c in range(3):ax.text(c,r,str(arr[r,c]),ha='center',va='center',fontsize=9,color='white' if arr[r,c]>185 else navy)
ax.set_xticks([0,1,2],['Input 1','Input 2','Input 3']);ax.set_yticks(range(9),[f'Output {i}' for i in range(1,10)])
ax.tick_params(length=0);ax.set_title('Barrier pairs in each\nnormalization scenario',loc='left',pad=14,fontsize=10,color=navy)
ax.text(-.30,1.10,'c',transform=ax.transAxes,fontsize=14,weight='bold',color=navy)
for spine in ax.spines.values():spine.set_visible(False)
fig.text(.645,.065,'Point summary: 73 pairs. Envelope: 0 pairs.\nAll counts use the same 2,556 endpoint pairs.\nShared preparations; dependent scenarios.',fontsize=9,linespacing=1.4,color=navy)
fig.text(.085,.016,'Empirical ranges are not confidence intervals. No independent validation or certified barrier is established.',fontsize=9,color=navy)
pdf=args.out/'Figure_FAS_Ordinal.pdf';png=args.out/'FAS_Ordinal_Evidence.png'
fig.savefig(pdf,metadata={'Title':'Measured FAS order constraints and conditional accessibility','Author':'Vivaan Patni','Creator':'','Producer':'','CreationDate':None,'ModDate':None})
fig.savefig(png,dpi=220,metadata={'Author':'Vivaan Patni'})
plt.close(fig)
(args.out/'figure_provenance.json').write_text(json.dumps({'source_sha256':hashlib.sha256(raw).hexdigest(),'figure_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'all_vertices':72,'all_primary_pairs':2556,'all_dependent_scenarios':27,'selection':'No outcome subset. Five comparisons are the first frozen solver certificate; no minimality or optimized robustness claimed.'},indent=2)+'\n',encoding='utf-8')
print(json.dumps({'pdf':str(pdf),'png':str(png)}))
