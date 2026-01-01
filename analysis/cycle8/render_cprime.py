"""Plot every declared B1 group and sensitivity from verified exact results."""
import argparse, gzip, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

p=argparse.ArgumentParser()
p.add_argument('--result',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
raw=a.result.read_bytes()
if a.result.suffix=='.gz':raw=gzip.decompress(raw)
expected='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'
if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('Result differs from independently verified primary')
d=json.loads(raw)
if d['status']!='COMPLETE_DESCRIPTIVE_AUDIT':raise ValueError('Incomplete source audit')
a.out.mkdir(parents=True,exist_ok=True)
primary=d['sensitivities'][0]
ids=sorted({r['intended'] for r in primary['group_cell_results']})
if len(ids)!=30:raise ValueError('Changed cohort')
by={(r['intended'],r['cell']):r for r in primary['group_cell_results']}
navy='#18354B';teal='#007F86';orange='#C07A20';grey='#96A5AF';light='#E7EEF1'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.linewidth':.7})
fig=plt.figure(figsize=(8.5,8.3),facecolor='white')
fig.text(.055,.971,'cPRIME: marker sensitivity and context limits',fontsize=14,weight='bold',color=navy)
fig.text(.055,.941,'30 intended labels  |  66 local states  |  396 paired measurements',fontsize=10,color=navy)
labels=[]
for s in ids:
 _,pos,ref,alt=s.split('-');labels.append(f'{pos}  {ref}>{alt}')
axes=[]
for ci,cell in enumerate(['K562','KBM7']):
 ax=fig.add_axes([.225+ci*.39,.345,.345,.535]);axes.append(ax)
 for i,s in enumerate(ids):
  r=by[(s,cell)]
  lo,hi,pt=[math.log2(float(F(r[k]))) for k in ('minimum_fold','maximum_fold','point_fold')]
  if not -1e-12<=lo<=pt+1e-12<=hi+2e-12:raise ValueError('Bounds or median outside expected order')
  if i%2==0:ax.axhspan(i-.5,i+.5,color='#F4F7F8',lw=0,zorder=0)
  ax.plot([lo,hi],[i,i],color=grey,lw=1.4,zorder=2)
  ax.plot([lo,lo],[i-.18,i+.18],color=orange if F(r['minimum_fold'])>1 else grey,lw=1.7)
  ax.plot([hi,hi],[i-.18,i+.18],color=grey,lw=1.2)
  ax.scatter([pt],[i],s=11,color=teal,zorder=4)
 ax.axvline(1,color=orange,lw=.8,ls=':',zorder=1)
 ax.set_xlim(-.12,8.4);ax.set_ylim(29.6,-.6)
 ax.set_xticks([0,2,4,6,8],['1','4','16','64','256'],fontsize=10)
 ax.set_yticks(range(30),labels if ci==0 else [],fontsize=10)
 ax.tick_params(axis='y',length=0,pad=8)
 ax.spines['left'].set_visible(False)
 ax.grid(axis='x',color=light,lw=.6,zorder=0)
 ax.set_title(cell,fontsize=11,weight='bold',color=navy,pad=9)
 ax.set_xlabel('Marker dispersion (fold; log scale)',labelpad=5,fontsize=10)
fig.text(.035,.898,'a',fontsize=14,weight='bold',color=navy)
fig.text(.072,.898,'Intended label (chr17)',fontsize=10,color=navy,va='bottom')
handles=[Line2D([0],[0],color=grey,lw=1.4,marker='|',markersize=7,label='Exact allowed range'),
 Line2D([0],[0],color=teal,lw=0,marker='o',markersize=4,label='Median-based dispersion'),
 Line2D([0],[0],color=orange,lw=.8,ls=':',label='Twofold display threshold')]
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.52,.278),ncol=3,frameon=False,fontsize=10,
 columnspacing=1.1,handletextpad=.6)

ax=fig.add_axes([.06,.085,.43,.135]);ax.axis('off')
ax.text(-.01,1.16,'b',fontsize=14,weight='bold',color=navy,transform=ax.transAxes)
ax.text(.085,1.17,'All 60 group/cell comparisons',fontsize=10,weight='bold',color=navy,transform=ax.transAxes)
table=[]
for sc in d['sensitivities']:
 s=sc['summary']
 # Count independently from saved exact group records, then check the display totals.
 records=sc['group_cell_results']
 counts=[sum(F(r['minimum_fold'])>1 for r in records),sum(F(r['minimum_fold'])>=2 for r in records),
         sum(F(r['maximum_fold'])<=2 for r in records)]
 table.append([f'{float(F(sc["gamma"])):g}',*map(str,counts)])
t=ax.table(cellText=table,colLabels=['Widening','Min > 1','Min ≥ 2','Max ≤ 2'],
 cellLoc='center',bbox=[0,0,1,1],colWidths=[.26,.245,.245,.25])
t.auto_set_font_size(False);t.set_fontsize(10)
for (row,col),cell in t.get_celld().items():
 cell.set_linewidth(.4);cell.set_edgecolor('white')
 cell.set_facecolor(light if row==0 else ('#F4F7F8' if row%2 else 'white'))
 cell.get_text().set_color(navy)
 if row==0:cell.get_text().set_weight('bold')

ax=fig.add_axes([.615,.085,.335,.135])
exons=[17,25,27,39,42]
for j,sc in enumerate(d['sensitivities']):
 records={int(r['exon']):r for r in sc['exon_shared_order_results']}
 for i,ex in enumerate(exons):
  excluded=ex in sc['summary']['shared_order_infeasible_exons']
  ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=orange if excluded else light,edgecolor='white',lw=2))
  ax.text(j,i,'E' if excluded else 'N',ha='center',va='center',fontsize=10,color='white' if excluded else navy)
ax.set_xlim(-.5,3.5);ax.set_ylim(4.5,-.5)
ax.set_xticks(range(4),['1','1.25','1.5','2'],fontsize=10)
ax.set_yticks(range(5),[str(x) for x in exons],fontsize=10)
ax.tick_params(length=0);ax.set_ylabel('Exon',fontsize=10,labelpad=6)
ax.set_xlabel('Widening factor',fontsize=10,labelpad=3)
for sp in ax.spines.values():sp.set_visible(False)
ax.text(-.27,1.16,'c',fontsize=14,weight='bold',color=navy,transform=ax.transAxes)
ax.text(-.14,1.17,'Common order across cells',fontsize=10,weight='bold',color=navy,transform=ax.transAxes)
fig.text(.058,.014,'Preparation ranges are descriptive, not confidence intervals.',fontsize=10,color=navy)
fig.text(.615,.014,'E: excluded     N: not excluded',fontsize=10,color=navy)
fig.savefig(a.out/'Figure_cPRIME.pdf',metadata={'Title':'Marker sensitivity and context limits','Author':'Vivaan Patni','Creator':'','Producer':'','CreationDate':None,'ModDate':None})
fig.savefig(a.out/'Figure_cPRIME.png',dpi=200,metadata={'Author':'Vivaan Patni'})
receipt={'status':'generated_awaiting_visual_review','primary_sha256':expected,'groups':30,'displayed_ranges':60,
 'gamma_count':4,'context_tests':20,'minimum_native_font_points':10,
 'figure_inches':[8.5,8.3], 'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(a.out.glob('Figure_cPRIME.*'))}}
(a.out/'Figure_Provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
