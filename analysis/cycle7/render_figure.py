"""Display all frozen successor diagnostics, without fitting or selecting outcomes."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--robustness',type=Path,required=True)
p.add_argument('--geometry',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
def read(path,expected):
    raw=path.read_bytes()
    if path.suffix=='.gz':raw=gzip.decompress(raw)
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('Unverified input')
    return json.loads(raw)
r=read(a.robustness,'861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa')
g=read(a.geometry,'0c2ad1a6a602f80e44ff05eaf3321522d9eaf0940382c599309a04d33f4b383e')
navy='#18354B';teal='#007F86';orange='#D18A24';red='#B34D58';grey='#C4CDD3'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,
    'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.linewidth':.7})
fig=plt.figure(figsize=(8.5,8.3),facecolor='white')
fig.text(.09,.965,'FAS: limits of model exclusion and route certainty',fontsize=14,weight='bold',color=navy)
fig.text(.09,.936,'Exploratory follow-up of the complete, frozen 72-state analysis',fontsize=10,color=navy)
for x,value,label in [(.09,'1.029716','smallest model-rescuing factor'),(.43,'0.016724','smallest additive expansion'),(.73,'0 / 2,556','all-scenario obstructions')]:
    fig.text(x,.889,value,fontsize=17,weight='bold',color=teal)
    fig.text(x,.863,label,fontsize=9,color=navy)
gs=fig.add_gridspec(2,2,left=.115,right=.965,bottom=.25,top=.775,width_ratios=[1.4,1],height_ratios=[1,1.14],wspace=.40,hspace=.70)
ax=fig.add_subplot(gs[0,0]);full=float(F(r['full']['multiplicative_squared']['gamma_enclosure']['lower']))-1
vals=[100*(float(F(z['results']['multiplicative_squared']['gamma_enclosure']['lower']))-1) for z in r['all_single_preparation_deletions']]
labels=[f"{z['kind'].capitalize()} {z['omitted_zero_based']+1}" for z in r['all_single_preparation_deletions']]
y=np.arange(12)
ax.axvline(100*full,color=navy,lw=1,ls=':',label='All preparations')
ax.hlines(y,100*full,vals,color=grey,lw=1)
ax.scatter(vals,y,s=23,color=teal,zorder=3)
for v,k in zip(vals,y):ax.annotate(f'{v:.2f}',(v,k),xytext=(5,0),textcoords='offset points',va='center',fontsize=9,color=navy)
ax.set_yticks(y,labels,fontsize=9);ax.invert_yaxis();ax.set_xlim(0,17.4)
ax.set_xticks([0,5,10,15]);ax.set_xlabel('Expansion above one (%)')
ax.set_title('One preparation omitted at a time',loc='left',pad=11,fontsize=10.5,color=navy)
ax.text(-.16,1.10,'a',transform=ax.transAxes,fontsize=14,weight='bold',color=navy)
ax.grid(axis='x',color='#E5E9EC',lw=.5);ax.legend(frameon=False,fontsize=9,loc='lower right')

ax=fig.add_subplot(gs[0,1])
for z in r['all_single_vertex_deletions']:
    v=z['omitted_vertex'];c=navy if v==46 else orange if z['certificate']['status']=='feasible' else grey
    ax.scatter(v%12,5-v//12,s=40,color=c,marker='s')
ax.set(xlim=(-.6,11.6),ylim=(-5,6));ax.axis('off')
ax.set_title('All 72 observation omissions',loc='left',pad=11,fontsize=10.5,color=navy)
ax.text(-.09,1.10,'b',transform=ax.transAxes,fontsize=14,weight='bold',color=navy)
ax.legend(handles=[Patch(facecolor=grey,label='69: exclusion persists'),Patch(facecolor=orange,label='13, 42: model becomes feasible'),Patch(facecolor=navy,label='46: fixed reference identity')],
    loc='lower left',bbox_to_anchor=(-.05,-.09),frameon=False,fontsize=9,labelspacing=.7)

ax=fig.add_subplot(gs[1,:])
exc=[z for z in g['pairs'] if F(z['C_fix'])<F(z['E'])]
if len(exc)!=81:raise ValueError('Unexpected pair inventory')
for k,z in enumerate(exc):
    e=F(z['E']);q=[0,float(F(z['C_fix'])/e),float(F(z['C_ad'])/e),1]
    if q[1]<.8:raise ValueError('Display zoom hides a non-fixed interval')
    for start,end,c in zip(q[:-1],q[1:],[teal,orange,red]):
        ax.barh(k+1,end-start,left=start,height=1,color=c,edgecolor='none')
ax.set_xlim(.8,1);ax.set_ylim(81.5,.5);ax.set_yticks([1,21,41,61,81])
ax.set_xticks([.8,.85,.9,.95,1]);ax.set_xlabel(r'Common threshold / endpoint ceiling ($\tau/E_{s,t}$), upper 20% shown')
ax.set_ylabel('Endpoint-pair rank\n(fixed genotype order)')
ax.set_title('Every pair with a non-fixed route interval (81 of 2,556)',loc='left',pad=12,fontsize=11,color=navy)
ax.text(-.075,1.10,'c',transform=ax.transAxes,fontsize=14,weight='bold',color=navy)
ax.legend(handles=[Patch(facecolor=teal,label='One route in all scenarios'),Patch(facecolor=orange,label='A route in each; no common route'),Patch(facecolor=red,label='Some scenarios connected')],
    loc='upper left',bbox_to_anchor=(-.01,-.38),frameon=False,fontsize=9,ncol=1)
fig.text(.54,.095,'2,475 pairs have a fixed route throughout.\n72 have adaptive-only intervals; 22 have mixed\nintervals (13 pairs have both). No pair is blocked\nin all scenarios at a common eligible threshold.',fontsize=9,color=navy,linespacing=1.4)
fig.text(.09,.035,'Empirical envelopes are sensitivity sets, not confidence intervals. Preparation omissions narrow them.\nThe factor rescues one box-compatible field; it is not measured recovery bias or a correction of all 27 fields.',fontsize=9,color=navy,linespacing=1.4)
pdf=a.out/'Figure_FAS_Robustness.pdf';png=a.out/'FAS_Robustness_and_Routes.png'
fig.savefig(pdf,metadata={'Title':'FAS robustness and joint-scenario route sensitivity','Author':'Vivaan Patni','Creator':'','Producer':'','CreationDate':None,'ModDate':None})
fig.savefig(png,dpi=220,metadata={'Author':'Vivaan Patni'})
plt.close(fig)
(a.out/'figure_provenance.json').write_text(json.dumps({'figure':pdf.name,'input_raw_sha256':{'robustness':'861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa','scenario_geometry':'0c2ad1a6a602f80e44ff05eaf3321522d9eaf0940382c599309a04d33f4b383e'},'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'png_sha256':hashlib.sha256(png.read_bytes()).hexdigest(),'status':'generated_pending_visual_review'},indent=2)+'\n')
print(pdf)
