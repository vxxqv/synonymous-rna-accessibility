"""Fresh evidence-led figures; figure assets must remain outside this archive."""
from pathlib import Path
import json,csv,argparse,math,sys,itertools,textwrap,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
ROOT=Path(__file__).resolve().parents[1]
INK='#20272B';TEAL='#176B68';BLUE='#42699B';ORANGE='#B85E31';PURPLE='#78628C';GREY='#8A939B';LIGHT='#DCE1E4'
CONDITIONS=('incell','cellfree');LABELS={'incell':'In cell','cellfree':'Cell free'};COLORS={'incell':TEAL,'cellfree':BLUE}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'xtick.labelsize':8,
 'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.2,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','text.color':INK,'axes.labelcolor':INK,'legend.frameon':False})
SEEN=set()
def load(path):
    path=Path(path);SEEN.add(path.as_posix());return json.loads((ROOT/path).read_text())
def family(name):return [load(p.relative_to(ROOT)) for p in sorted((ROOT/'results'/name).glob('*.json')) if p.name!='summary.json']
def panel(ax,letter,title):
    ax.text(-.15,1.08,letter,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
    ax.set_title(textwrap.fill(title,30),loc='left',pad=9,fontweight='bold',fontsize=8.5)
def grid(rows=3,cols=2,height=153):
    fig,axes=plt.subplots(rows,cols,figsize=(174/25.4,height/25.4),squeeze=False)
    fig.subplots_adjust(left=.17,right=.975,top=.91,bottom=.095,wspace=.65,hspace=1.05)
    return fig,axes
def ecdf(ax,values,**kwargs):
    v=np.sort(values);ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',**kwargs)
def interval(ax,row,y,color=TEAL,marker='o'):
    ax.errorbar(row['estimate'],y,xerr=[[row['estimate']-row['lower']],[row['upper']-row['estimate']]],fmt=marker,color=color,ms=4,capsize=2,lw=.9)
def save(fig,name,out,data):
    fig.canvas.draw()
    for ax in fig.axes:
        for axis,bounds in ((ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())):
            lower,upper=sorted(bounds)
            for tick in axis.get_major_ticks():
                if not lower-1e-10<=tick.get_loc()<=upper+1e-10:tick.label1.set_visible(False);tick.label2.set_visible(False)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();bbox=fig.bbox
    clipped=[]
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0<bbox.x0-1 or b.y0<bbox.y0-1 or b.x1>bbox.x1+1 or b.y1>bbox.y1+1:clipped.append(t.get_text())
    if clipped:raise ValueError(name+': clipped text '+str(clipped))
    data['export']={'width_mm':174,'height_mm':float(fig.get_figheight()*25.4),'clipped_text_preflight':clipped,'ordinary_font_pt':8.5}
    data['provenance']={'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs':[{'path':p,'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in sorted(SEEN)],'matplotlib':matplotlib.__version__,'numpy':np.__version__}
    (ROOT/'results/figure_data'/f'{name}.json').write_text(json.dumps(data,indent=2,allow_nan=False))
    for ext in ('pdf','svg','eps'):fig.savefig(out/f'{name}.{ext}',dpi=600)
    fig.savefig(out/f'{name}.png',dpi=600);plt.close(fig)
    print(name,'text outside canvas',clipped,flush=True)

def figure1(out):
    filt=family('expansion/filtration');smooth=family('expansion/smoothness');add=family('expansion/additive');inf=load('results/expansion/inference.json')
    fig,axes=grid(height=158);data={'claim':'Most contrast against shuffled scores is generic organization; additive surrogates are similarly connected','population':'404 complete blocks,185 transcript contexts; discovery cohort','panels':{}}
    ax=axes[0,0];panel(ax,'a','All-threshold fragmentation')
    use=[r for r in filt if r['condition']=='incell'];suse=[r for r in smooth if r['condition']=='incell'];ause=[r for r in add if r['condition']=='incell']
    series=[('Folding score',[r['observed']['integrated_excess'] for r in use],TEAL,'-'),('Score shuffle',[r['null']['integrated_excess']['mean'] for r in use],GREY,':'),('Smooth reference',[r['heat']['integrated_excess']['mean'] for r in suse],BLUE,'--'),('Additive codon',[r['approximations']['codon']['rank_restored']['integrated_excess']['mean'] for r in ause],PURPLE,'-.')]
    for label,values,col,ls in series:ecdf(ax,values,color=col,ls=ls,label=label)
    ax.set(xlabel=r'Integrated excess  $I_E$',ylabel='Cumulative fraction',xlim=(0,.24),ylim=(0,1.03));ax.legend(loc='lower right',handlelength=1.6,fontsize=7.2,labelspacing=.2,borderaxespad=.3)
    data['panels']['a']={k:v for k,v,_,_ in series}
    ax=axes[0,1];panel(ax,'b','Observed versus smooth reference')
    for c,mark in zip(CONDITIONS,('o','x')):
        rr=[r for r in smooth if r['condition']==c];x=[r['observed']['integrated_excess'] for r in rr];y=[r['heat']['integrated_excess']['mean'] for r in rr]
        ax.scatter(x,y,c=COLORS[c],s=11,marker=mark,lw=.5,label=LABELS[c],rasterized=True);data['panels'][f'b_{c}']={'observed':x,'reference':y}
    ax.plot([0,.32],[0,.32],color=GREY,ls=':',lw=.8);ax.set(xlabel=r'Folding-score $I_E$',ylabel='Smooth-reference\nmean $I_E$',xlim=(-.005,.32),ylim=(-.005,.32));ax.legend(loc='upper left')
    ax=axes[1,0];panel(ax,'c','Residual contrast after organization')
    keys=[('smoothness','integrated_excess','Heat-rank'),('additive','nucleotide','Additive base'),('additive','codon','Additive codon')];rows=[]
    for j,(src,key,label) in enumerate(keys):
        for k,c in enumerate(CONDITIONS):
            row=inf[f'{src}/{c}/{key}/null_minus_observed'];interval(ax,row,j+(k-.5)*.18,COLORS[c],('o','s')[k]);rows.append({'reference':label,'condition':c,**row})
    ax.axvline(0,color=GREY,lw=.8,ls=':');ax.set(yticks=range(3),yticklabels=[k[2] for k in keys],xlabel=r'Reference − folding-score mean $I_E$',ylim=(-.6,2.7));ax.invert_yaxis();data['panels']['c']=rows
    ax=axes[1,1];panel(ax,'d','Degree matching at cutoff 0.5')
    rows=[]
    for j,name in enumerate(('uniform','degree')):
        for k,c in enumerate(CONDITIONS):
            row=inf[f'conditional_null/{c}/{name}/null_minus_observed'];interval(ax,row,j+(k-.5)*.18,COLORS[c],('o','s')[k]);rows.append({'reference':name,'condition':c,**row})
    ax.axvline(0,color=GREY,lw=.8,ls=':');ax.set(yticks=[0,1],yticklabels=['Size matched','Degree + size'],xlabel=r'Reference − observed $E_{0.5}$',xlim=(0,.14),ylim=(-.6,1.6));ax.invert_yaxis();data['panels']['d']=rows
    ax=axes[2,0];panel(ax,'e','Exact roughness-conditioned control')
    rows=[]
    for k,c in enumerate(CONDITIONS):
        for j,t in enumerate(('0.05','0.1','0.2')):
            row=inf[f'smoothness/{c}/exact_conditional/{t}'];interval(ax,row,j+(k-.5)*.18,COLORS[c],('o','s')[k]);rows.append({'tolerance':t,'condition':c,**row})
    ax.axvline(0,color=GREY,lw=.8,ls=':');ax.set(yticks=range(3),yticklabels=['±5%','±10%','±20%'],xlabel=r'Reference − observed mean $I_E$',ylim=(-.6,2.7));ax.invert_yaxis();data['panels']['e']=rows
    ax=axes[2,1];panel(ax,'f','Additive approximation is descriptive')
    for c,mark in zip(CONDITIONS,('o','x')):
        rr=[r for r in add if r['condition']==c];x=[r['approximations']['codon']['in_sample_R_squared'] for r in rr];y=[r['observed']['integrated_excess']-r['approximations']['codon']['rank_restored']['integrated_excess']['mean'] for r in rr]
        ax.scatter(x,y,s=11,c=COLORS[c],marker=mark,lw=.5,rasterized=True);data['panels'][f'f_{c}']={'in_sample_R_squared':x,'observed_minus_additive_I_E':y}
    ax.axhline(0,color=GREY,lw=.8,ls=':');ax.set(xlabel=r'Codon-additive in-sample $R^2$',ylabel=r'Folding − additive $I_E$',xlim=(-.03,1.03));save(fig,'Figure_1',out,data)

def figure3(out):
    exact=load('analysis/exact/exact_peptide_automata.json');cost=load('analysis/exact/exact_bridge_threshold_automata.json');fig,axes=grid(2,2,129);axes=axes.flat;data={'claim':'Codon register phase determines exact finite-language scarcity, not a physical phase transition','panels':{}}
    styles=('-','--',':');cols=(TEAL,BLUE,ORANGE);lam=[25.616347998891122,17.96947880175194,28.15737830833643];mu=[9.,9.518542991203905,11.515466277638433]
    for d,col,ls in zip((1,2,3),cols,styles):
        b=np.array([int(r['both_nonempty']) for r in exact[str(d)]['exact_counts'][:40]],dtype=object);accepted=np.array([int(r['delta_le1_peptide_pairs']) for r in cost[str(d)]['exact_counts'][:40]],dtype=object);m=np.arange(1,len(b)+1);p=np.array([int(n)/84**int(k) for n,k in zip(b,m)]);q=np.array([int(u)/int(v) for u,v in zip(accepted,b)])
        axes[0].semilogy(m,p,color=col,ls=ls,label=f'Shift {d}');axes[1].semilogy(m,q,color=col,ls=ls)
        local=np.diff(np.log(q));axes[2].plot(m[1:],local,color=col,ls=ls);axes[2].axhline(math.log(mu[d-1]/lam[d-1]),color=col,ls=':',lw=.65)
        data['panels'][str(d)]={'codons_per_arm':m.tolist(),'both_feasible_count':[str(v) for v in b],'near_transition_count':[str(v) for v in accepted],'conditional_feasibility':p.tolist(),'conditional_near_transition':q.tolist(),'lambda':lam[d-1],'mu':mu[d-1]}
    for ax,l,t in zip(axes,'abcd',('Shifted-target feasibility','Zero-or-one-mutation transitions','Finite-length logarithmic rate','Certified exponential roots')):panel(ax,l,t)
    axes[0].set(xlabel='Codons per arm  m',ylabel='Conditional probability (log)');axes[0].legend(loc='lower left')
    axes[1].set(xlabel='Codons per arm  m',ylabel='Conditional probability (log)')
    axes[2].set(xlabel='Codons per arm  m',ylabel=r'$\log(q_m/q_{m-1})$')
    for j,(l,u,col) in enumerate(zip(lam,mu,cols)):
        axes[3].plot([u,l],[j,j],color=LIGHT,lw=2);axes[3].plot(l,j,'o',color=col,ms=5);axes[3].plot(u,j,'s',mfc='white',mec=col,ms=5)
        axes[3].text(l+.5,j,f'{l:.3f}',va='center',fontsize=8)
    axes[3].set(yticks=range(3),yticklabels=['Shift 1','Shift 2','Shift 3'],xlabel='Integer-matrix growth root',xlim=(0,34));axes[3].invert_yaxis();axes[3].legend([Line2D([],[],marker='o',ls='',color=INK),Line2D([],[],marker='s',ls='',mfc='white',color=INK)],[r'Both targets $\lambda_d$',r'Near transition $\mu_d$'],loc='lower left',bbox_to_anchor=(-.05,1.22),ncol=2,fontsize=7.2,handlelength=.8,columnspacing=.8)
    save(fig,'Figure_3',out,data)

def figure2(out):
    audit=load('results/expansion/verification_release/arbitrary_pair_graph_results.json');source=load('results/expansion/verification_release/general_graph_figure_sources.json')
    fig=plt.figure(figsize=(174/25.4,157/25.4));gs=fig.add_gridspec(3,2,left=.105,right=.97,top=.91,bottom=.105,wspace=.7,hspace=1.05)
    data={'claim':'Complete standard-code component signatures coexist with variable path multiplicity and code-induced completion weights','universe':'all simple pairing graphs on one or two codons; abstract compatibility, not physical multi-partner structures','panels':{}}
    ax=fig.add_subplot(gs[:2,0]);panel(ax,'a','Components across 400 peptides');aa=sorted(set(r['peptide'][0] for r in audit['peptide_component_summaries'] if r['codons']==2));mapping={r['peptide']:r for r in audit['peptide_component_summaries'] if r['codons']==2};matrix=np.array([[mapping[a+b]['maximum_components'] for b in aa] for a in aa]);im=ax.imshow(matrix,vmin=1,vmax=4,cmap=LinearSegmentedColormap.from_list('classes',['#E8ECEE','#94BBB6','#176B68']),interpolation='nearest')
    ax.set(xticks=range(20),xticklabels=aa,yticks=range(20),yticklabels=aa,xlabel='Second amino acid',ylabel='First amino acid');ax.tick_params(length=0)
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.05,pad=.21,ticks=[1,2,3,4]);cb.set_label('Maximum over all compatible pair graphs',fontsize=8);data['panels']['a']={'amino_acid_order':aa,'maximum_component_matrix':matrix.tolist(),'scope':'all400 ordered two-amino-acid peptides, every simpleH on six nucleotide positions'}
    ax=fig.add_subplot(gs[0,1]);panel(ax,'b','Exact distance distribution');counts=[r for r in audit['distance_counts'] if r['status']=='connected'];x=[r['hamming_distance'] for r in counts];y=[r['endpoint_pairs'] for r in counts];ax.bar(x,np.array(y)/1000,color=TEAL,width=.6)
    ax.set(xlabel='Graph distance = Hamming distance',ylabel='Endpoint pairs (thousands)',xticks=x,ylim=(0,510));ax.text(.5,.96,'719,338 pairs; zero violations',transform=ax.transAxes,ha='center',va='top',fontsize=8);data['panels']['b']=counts
    ax=fig.add_subplot(gs[1,1]);panel(ax,'c','Permitted mutation-order fraction');fractions=sorted(audit['mutation_order_fraction_counts'],key=lambda v:float(__import__('fractions').Fraction(v)),reverse=True);counts=[audit['mutation_order_fraction_counts'][v] for v in fractions];ax.plot(counts,range(len(counts)),'o',c=TEAL,ms=4);ax.set_xscale('log');ax.set(yticks=range(len(counts)),yticklabels=fractions,xlabel='Endpoint pairs (log scale)',ylabel=r'$\mathsf{P}(x,y)/k!$',ylim=(-.6,len(counts)-.4));ax.invert_yaxis();data['panels']['c']={'order_fractions':fractions,'endpoint_counts':counts}
    ax=fig.add_subplot(gs[2,0]);panel(ax,'d','A code-induced volume prediction');star=source['star_family'][:16];r=[s['partners'] for s in star];p=[s['C_first_probability_numerator']/s['C_first_probability_denominator'] for s in star];ax.semilogy(r,p,'o-',c=TEAL,ms=3);ax.set(xlabel='Required partner constraints  r',ylabel='Leu-first C probability',xticks=[1,5,10,15]);data['panels']['d']=star
    ax=fig.add_subplot(gs[2,1]);panel(ax,'e','Branching changes route counts');star=source['star_family'][:20];k=[s['changed_events'] for s in star];counts=[int(s['shortest_paths']) for s in star];ax.semilogy(k,counts,'o-',c=BLUE,ms=3);ax.set(xlabel='Differing nucleotide positions',ylabel='Exact shortest routes (log)');ax.text(.04,.91,r'$\mathsf{P}=r!$ for $\mathrm{LK}^{r}$',transform=ax.transAxes,fontsize=9);data['panels']['e']=star
    data['audit_counts']=audit['counts'];save(fig,'Figure_2',out,data)

def field_grid(row,thresholds):
    curve=row['observed']['curve'];levels=np.array([r['threshold'] for r in curve]);scores=np.array([r['excess'] for r in curve]);out=[]
    for t in thresholds:
        candidates=np.flatnonzero(levels>=t);out.append(float(scores[candidates[-1]]) if len(candidates) else 0.)
    return out

def figure4(out):
    rows=family('expansion/filtration');fig,axes=grid(height=159);fig.subplots_adjust(right=.89);data={'claim':'A cutoff can miss barriers, while most observed barriers have short retention-threshold lifetimes','panels':{}};thresholds=np.linspace(0,1,101)
    by={(r['id'],r['condition']):r for r in rows};ids=sorted(set(r['id'] for r in rows),key=lambda i:-max(by[i,c]['observed']['longest_lifetime'] for c in CONDITIONS));cmap=LinearSegmentedColormap.from_list('excess',['#FFFFFF','#6EAAA2','#134D4B'])
    affected=[i for i in ids if any(by[i,c]['observed']['positive_bar_count']>0 for c in CONDITIONS)]
    for j,c in enumerate(CONDITIONS):
        ax=axes[0,j];n=sum(by[i,c]['observed']['positive_bar_count']>0 for i in ids);panel(ax,chr(97+j),f'{LABELS[c]}: {n}/404 with a bar');all_values=np.array([field_grid(by[i,c],thresholds) for i in ids]);values=np.array([field_grid(by[i,c],thresholds) for i in affected]);im=ax.imshow(values,aspect='auto',origin='upper',extent=(0,1,len(affected),0),cmap=cmap,vmin=0,vmax=1,interpolation='nearest',rasterized=True)
        ax.axvline(.5,color=INK,lw=.65,ls=':');ax.set(xlabel='Retention threshold',ylabel='Shared affected-block rank',yticks=[1,len(affected)//2,len(affected)]);data['panels'][c]={'all_ordered_ids':ids,'displayed_affected_ids':affected,'threshold':thresholds.tolist(),'all404_E_F':all_values.tolist(),'displayed_E_F':values.tolist(),'omitted_zero_bar_blocks':404-len(affected)}
    cb=fig.colorbar(im,ax=list(axes[0]),shrink=.7,pad=.015,fraction=.025);cb.set_label(r'Extra fragmentation $E_F$',fontsize=8)
    ax=axes[1,0];panel(ax,'c','Persistence across all thresholds')
    for c,ls in zip(CONDITIONS,('-','--')):
        use=[r for r in rows if r['condition']==c];values=[r['observed']['longest_lifetime'] for r in use];ecdf(ax,values,color=COLORS[c],ls=ls,label=LABELS[c]);data['panels'][f'c_{c}']=values
    ax.set(xlabel='Longest positive lifetime',ylabel='Cumulative fraction',ylim=(0,1.03));ax.legend(loc='lower right')
    ax=axes[1,1];panel(ax,'d','Barrier depth and total effect')
    for c,marker in zip(CONDITIONS,('o','x')):
        use=[r for r in rows if r['condition']==c];x=[r['observed']['longest_lifetime'] for r in use];y=[r['observed']['integrated_excess'] for r in use];ax.scatter(x,y,s=12,c=COLORS[c],marker=marker,lw=.5,rasterized=True);data['panels'][f'd_{c}']={'lifetime':x,'integral':y}
    ax.set(xlabel='Longest positive lifetime',ylabel=r'Integrated excess $I_E$')
    ax=axes[2,0];panel(ax,'e','Largest observed persistence intervals');bars=sorted([{'id':r['id'],'gene':r['gene'],'condition':r['condition'],**b} for r in rows for b in r['observed']['bars']],key=lambda r:-r['lifetime'])[:10]
    for i,b in enumerate(bars):ax.plot([b['death'],b['birth']],[i,i],color=COLORS[b['condition']],lw=2.2);ax.plot(b['death'],i,'|',c=COLORS[b['condition']],ms=5)
    ax.set(yticks=range(len(bars)),yticklabels=[b['gene']+(' IC' if b['condition']=='incell' else ' CF') for b in bars],xlabel='Retention threshold',ylim=(-.7,len(bars)-.3),xlim=(0,1));ax.invert_yaxis();data['panels']['e']={'selection':'top10 finite bars by lifetime, descriptive exception display', 'bars':bars}
    ax=axes[2,1];panel(ax,'f','Counterfactual local repair');cases=load('results/local_gate_restoration.json');groups=[sum(r['fully_explained_by_single_GU_addbacks'] for r in cases),sum(r['fully_restored_by_GU_and_coding'] and not r['fully_explained_by_single_GU_addbacks'] for r in cases),sum(not r['fully_restored_by_GU_and_coding'] for r in cases)]
    ax.barh(range(3),groups,color=[TEAL,BLUE,ORANGE],height=.55)
    for i,n in enumerate(groups):ax.text(n+.35,i,str(n),va='center')
    ax.set(yticks=range(3),yticklabels=['GU addbacks','+ coding gates','Other squares'],xlabel='Cases fully repaired (n = 25)',xlim=(0,19));ax.invert_yaxis();data['panels']['f']={'mutually_exclusive_repair_counts':groups,'counterfactual':'Excluded genotypes are restored; this violates the original filter and does not establish traversability in vivo'}
    save(fig,'Figure_4',out,data)

def figure5(out):
    rs=family('expansion/rnastructure');long=family('expansion/longer');sens=load('results/sensitivity_summary.json');data={'claim':'Robustness checks support a bounded model result and expose limits; none supplies prospective biological validation','panels':{}};fig,axes=grid(height=158)
    ax=axes[0,0];panel(ax,'a','Independent engine, same sequences')
    for c,marker in zip(CONDITIONS,('o','x')):
        rr=[r for r in rs if r['condition']==c];x=[r['vienna_filtration']['integrated_excess'] for r in rr];y=[r['filtration']['integrated_excess'] for r in rr];ax.scatter(x,y,s=21,marker=marker,c=COLORS[c],lw=.8,label=LABELS[c]);data['panels'][f'a_{c}']={'ViennaRNA':x,'RNAstructure':y}
    lim=max([r['filtration']['integrated_excess'] for r in rs]+[r['vienna_filtration']['integrated_excess'] for r in rs])+.001;ax.plot([0,lim],[0,lim],color=GREY,ls=':',lw=.8);ax.set(xlabel=r'ViennaRNA $I_E$',ylabel=r'RNAstructure $I_E$',xlim=(-.0002,lim),ylim=(-.0002,lim));ax.legend(loc='upper left')
    ax=axes[0,1];panel(ax,'b','Longer complete variable blocks');native={(r['id'],r['condition']):r for r in family('expansion/filtration')}
    for n,dx in ((7,-.15),(9,0),(11,.15)):
        for r in long:
            if r['codons']!=n:continue
            base=native[r['id'],r['condition']]['observed']['integrated_excess'];value=r['filtration']['integrated_excess'];ax.plot([5+dx,n+dx],[base,value],color=COLORS[r['condition']],lw=.7);ax.plot(n+dx,value,'o',ms=3,c=COLORS[r['condition']])
    ax.set(xlabel='Variable codons',ylabel=r'Integrated excess $I_E$',xticks=[5,7,9,11]);data['panels']['b']=[{k:r[k] for k in ('id','condition','codons','n_compatible','cartesian_size','filtration')} for r in long]
    ax=axes[1,0];panel(ax,'c','Whole-threshold context sensitivity')
    names=('context210','context300','Turner1999');context_rows=[r for r in load('results/expansion/context_filtration.json')['rows'] if r['selection_group']=='hash_selected' and r['sensitivity'] in ('context210','context300')];contexts={(r['id'],r['condition'],r['sensitivity']):r for r in context_rows}
    selected_ids=sorted(set(r['id'] for r in context_rows))
    for i in selected_ids:
        for c in CONDITIONS:
            values=[native[i,c]['observed']['integrated_excess']]+[contexts[i,c,k]['filtration']['integrated_excess'] for k in ('context210','context300')];ax.plot([150,210,300],values,color=COLORS[c],lw=.6);ax.plot([150,210,300],values,'o',ms=2.5,c=COLORS[c])
    ax.set(xlabel='Fixed-context length (nt)',ylabel=r'Integrated excess $I_E$',xticks=[150,210,300]);data['panels']['c']=context_rows
    ax=axes[1,1];panel(ax,'d','Targeted exceptions are less stable')
    for j,k in enumerate(names):
        r=sens['exceptions'][k];counts=[r['fragmentation_persists'],r['nonempty']-r['fragmentation_persists'],r['empty']];left=0
        for value,color in zip(counts,(ORANGE,TEAL,LIGHT)):ax.barh(j,value,left=left,height=.56,color=color);left+=value
        ax.text(left+.4,j,f"{counts[0]}/{r['tested_original_cases']}",va='center')
    ax.set(yticks=range(3),yticklabels=['210 nt','300 nt','Turner 1999'],xlabel='Available original exception cases',xlim=(0,22));ax.invert_yaxis();data['panels']['d']=sens['exceptions']
    ax=axes[2,0];panel(ax,'e','Original pairing-sensitive source');pair=family('pairmap');rr=[r for r in pair if 'thresholds' in r];genes=['rpsM','rpsB'];vals=[]
    for j,gene in enumerate(genes):
        use=[r for r in rr if r['gene']==gene];v=[r['a_compatible'] for r in use];ax.scatter(np.array(v),np.full(len(v),j),s=24,c=[TEAL,BLUE][j],marker=['o','s'][j]);vals.append({'gene':gene,'condition_graphs':len(use),'a_compatible':v})
    ax.set(yticks=[0,1],yticklabels=['rpsM: 17 windows','rpsB: 3 windows'],xlabel='Compatibility connection probability',xlim=(.4,1.06),ylim=(-.6,1.6));ax.text(.03,.93,'Two loci; overlapping windows',transform=ax.transAxes,fontsize=7.2,ha='left',va='top');data['panels']['e']=vals
    ax=axes[2,1];panel(ax,'f','Evidence independence');sources=['Primary','PAIR-MaP','Contexts','Engine','Longer blocks'];matrix=[[0,0,0],[1,0,0],[0,0,0],[0,1,0],[0,0,0]]
    for j,row in enumerate(matrix):
        for k,value in enumerate(row):ax.plot(k,j,'o' if value else 'x',ms=5,color=TEAL if value else GREY,mfc=TEAL if value else 'white',mew=.8)
    ax.set(xticks=range(3),xticklabels=['New source','New engine','Prospective'],yticks=range(5),yticklabels=sources,xlim=(-.4,2.4),ylim=(5.7,-.7));ax.tick_params(axis='x',labelrotation=25);ax.text(1,5.25,'Prospective validation pending',ha='center',fontsize=8,color=ORANGE)
    ax.spines['left'].set_visible(False);ax.spines['bottom'].set_visible(False);ax.tick_params(length=0);data['panels']['f']={'sources':sources,'flags':matrix,'zero_marker':'x','one_marker':'filled circle'}
    fig.legend([Line2D([],[],lw=5,color=c) for c in (ORANGE,TEAL,LIGHT)],['Excess persists','Nonempty; no excess','Empty'],ncol=3,loc='lower center',bbox_to_anchor=(.55,.005),handlelength=.9,columnspacing=1.)
    save(fig,'Figure_5',out,data)

def supplement1(out):
    SEEN.add('results/native_metrics.csv')
    data={'claim':'Cohort composition and exact-space coverage constrain generalization','panels':{}};fig,axes=grid(2,2,126);metrics=list(csv.DictReader((ROOT/'results/native_metrics.csv').open()));blocks=load('data/native_windows.json');maps={r['transcript_id']:r for r in load('data/provenance/transcript_mapping.json')}
    ax=axes[0,0];panel(ax,'a','Complete local space coverage')
    for c,marker in zip(CONDITIONS,('o','x')):
        rr=[r for r in metrics if r['condition']==c];x=[int(r['cartesian_size']) for r in rr];y=[int(r['n_compatible']) for r in rr];ax.scatter(x,y,s=14,marker=marker,c=COLORS[c],lw=.5,rasterized=True);data['panels'][f'a_{c}']={'cartesian':x,'compatible':y}
    ax.set_xscale('log');ax.set_yscale('log');ax.set(xlabel='Cartesian sequences (log)',ylabel='Compatible sequences (log)')
    ax=axes[0,1];panel(ax,'b','Compatibility fragmentation')
    for c,col,ls in zip(('incell','cellfree','kasugamycin'),(TEAL,BLUE,ORANGE),('-','--',':')):
        x=[float(r['effective_components']) for r in metrics if r['condition']==c];ecdf(ax,x,color=col,ls=ls,label=LABELS.get(c,'Kasugamycin'));data['panels'][f'b_{c}']=x
    ax.set(xlabel='Effective components',ylabel='Cumulative fraction',ylim=(0,1.03));ax.legend(loc='lower right')
    ax=axes[1,0];panel(ax,'c','Genome-wide block locations');coords=[(maps[b['transcript_id']]['start']+b['transcript_pos0'] if b['strand']=='+' else maps[b['transcript_id']]['end']-b['transcript_pos0'])/1e6 for b in blocks];bins=np.linspace(0,4.64,31);ax.hist(coords,bins=bins,color=TEAL,edgecolor='white',lw=.5);ax.set(xlabel='U00096.2 coordinate (Mbp)',ylabel='Selected blocks',xlim=(0,4.64));data['panels']['c']={'coordinate_Mbp':coords,'bin_edges':bins.tolist()}
    ax=axes[1,1];panel(ax,'d','Pruning enables larger exact spaces');long=family('expansion/longer')
    for n,marker,col in ((7,'o',TEAL),(9,'s',BLUE),(11,'^',ORANGE)):
        rr=[r for r in long if r['codons']==n];ax.scatter([r['cartesian_size'] for r in rr],[r['n_compatible'] for r in rr],s=20,marker=marker,c=col,label=f'{n} codons');data['panels'][f'd_{n}']=[{k:r[k] for k in ('id','condition','cartesian_size','n_compatible')} for r in rr]
    ax.set_xscale('log');ax.set_yscale('log');ax.set(xlabel='Cartesian sequences (log)',ylabel='Listed compatible states (log)');ax.legend(loc='lower right');save(fig,'Figure_S1',out,data)

def supplement2(out):
    smooth=family('expansion/smoothness');add=family('expansion/additive');inf=load('results/expansion/inference.json');data={'claim':'Fit diagnostics, resource coverage and influence delimit the generic-score comparisons','panels':{}};fig,axes=grid(height=158)
    ax=axes[0,0];panel(ax,'a','Roughness calibration diagnostic')
    for c,marker in zip(CONDITIONS,('o','x')):
        rr=[r for r in smooth if r['condition']==c];x=[r['roughness_target'] for r in rr];y=[r['roughness_sample_mean'] for r in rr];ax.scatter(x,y,s=13,marker=marker,c=COLORS[c],lw=.5,rasterized=True);data['panels'][f'a_{c}']={'target':x,'reference_mean':y}
    ax.plot([0,40],[0,40],color=GREY,ls=':',lw=.7);ax.set_xscale('symlog',linthresh=.001);ax.set_yscale('symlog',linthresh=.001);ax.set(xlabel='Observed edge roughness (symlog)',ylabel='Reference roughness\n(symlog)')
    ax=axes[0,1];panel(ax,'b','Additive score variance explained')
    for c in CONDITIONS:
        for name,ls in (('nucleotide',':'),('codon','-')):
            x=[r['approximations'][name]['in_sample_R_squared'] for r in add if r['condition']==c];ecdf(ax,x,color=COLORS[c],ls=ls,label=LABELS[c]+(' base' if name=='nucleotide' else ' codon'));data['panels'][f'b_{c}_{name}']=x
    ax.set(xlabel=r'In-sample $R^2$',ylabel='Cumulative fraction',xlim=(0,1),ylim=(0,1.03));ax.legend(loc='upper left')
    ax=axes[1,0];panel(ax,'c','Additive fit degrees of freedom')
    for c,marker in zip(CONDITIONS,('o','x')):
        rr=[r for r in add if r['condition']==c];x=[r['n_compatible'] for r in rr];y=[r['approximations']['codon']['total_residual_df'] for r in rr];ax.scatter(x,y,s=13,marker=marker,c=COLORS[c],lw=.5,rasterized=True);data['panels'][f'c_{c}']={'compatible':x,'residual_df':y}
    ax.set(xlabel='Complete graph vertices',ylabel='Fit residual\ndegrees of freedom')
    ax=axes[1,1];panel(ax,'d','Leave-one-transcript influence')
    for c,ls in zip(CONDITIONS,('-','--')):
        rows=inf[f'smoothness/{c}/integrated_excess/leave_one_cluster_out']['details'];y=sorted(r['estimate'] for r in rows);ax.plot(range(1,len(y)+1),y,c=COLORS[c],ls=ls);data['panels'][f'd_{c}']=rows
    ax.axhline(0,c=GREY,ls=':',lw=.7);ax.set(xlabel='Exclusion estimate rank',ylabel='Heat − observed\nmean $I_E$')
    ax=axes[2,0];panel(ax,'e','Exact conditional-control coverage')
    for c,marker in zip(CONDITIONS,('o','x')):
        rr=[r for r in smooth if r['condition']==c];x=[r['n_compatible'] for r in rr];y=[math.log10(int(r['exact_conditional']['label_permutations'])) for r in rr];ax.scatter(x,y,s=13,marker=marker,c=COLORS[c],lw=.5,rasterized=True);data['panels'][f'e_{c}']={'vertices':x,'log10_label_permutations':y}
    ax.axhline(math.log10(40320),color=ORANGE,ls=':',lw=.9);ax.set_xscale('log');ax.set(xlabel='Graph vertices (log)',ylabel='Label permutations\n(log$_{10}$ count)')
    ax=axes[2,1];panel(ax,'f','Conditional controls can be invariant');rows=[]
    for j,c in enumerate(CONDITIONS):
        rr=[r for r in smooth if r['condition']==c and r['exact_conditional']['status']=='exact'];v=sum(r['exact_conditional']['tolerances']['0.05']['integrated_excess_varied'] for r in rr);ax.barh(j,v,color=COLORS[c],height=.5);ax.barh(j,len(rr)-v,left=v,color=LIGHT,height=.5);ax.text(len(rr)+1,j,f'{v}/{len(rr)}',va='center');rows.append({'condition':c,'varied':v,'exact_graphs':len(rr)})
    ax.set(yticks=[0,1],yticklabels=['In cell','Cell free'],xlabel='Exact graphs; labels: varied / all',xlim=(0,84));ax.invert_yaxis();data['panels']['f']=rows;save(fig,'Figure_S2',out,data)

def supplement3(out):
    cases=sorted(load('results/local_gate_restoration.json'),key=lambda r:(r['gene'],r['condition']));data={'claim':'All primary-cutoff exceptions have sufficient local addbacks, which violate the original filter','panels':{}};fig=plt.figure(figsize=(174/25.4,139/25.4));gs=fig.add_gridspec(2,2,left=.18,right=.97,top=.90,bottom=.13,wspace=.65,hspace=.9)
    ax=fig.add_subplot(gs[:,0]);panel(ax,'a','Complete 25-case repair audit');matrix=np.array([[r['fully_explained_by_single_GU_addbacks'],r['fully_restored_by_GU_and_coding'],r['fully_restored_by_all_single_addbacks']] for r in cases],int)
    for j,row in enumerate(matrix):
        for k,value in enumerate(row):ax.plot(k,j,'o' if value else 'x',ms=4.5,color=TEAL if value else GREY)
    abbreviations={'incell':'IC','cellfree':'CF','kasugamycin':'KAS'}
    ax.set(xticks=[0,1,2],xticklabels=['GU','+ codon','All local'],yticks=range(25),yticklabels=[r['gene']+' '+abbreviations[r['condition']] for r in cases],xlim=(-.5,2.5),ylim=(24.8,-.8));ax.tick_params(length=0);ax.spines['left'].set_visible(False);ax.spines['bottom'].set_visible(False);data['panels']['a']={'cases':cases,'full_restoration_flags':matrix.tolist()}
    ax=fig.add_subplot(gs[0,1]);panel(ax,'b','Original excess and endpoint size');x=[r['n_endpoints'] for r in cases];y=[r['original_excess'] for r in cases]
    for c,marker,col in zip(('incell','cellfree','kasugamycin'),('o','s','^'),(TEAL,BLUE,ORANGE)):
        rr=[r for r in cases if r['condition']==c];ax.scatter([r['n_endpoints'] for r in rr],[r['original_excess'] for r in rr],s=20,c=col,marker=marker,label=abbreviations[c])
    ax.set_xscale('log');ax.set(xlabel='Filtered endpoints (log)',ylabel=r'Original excess $E_{0.5}$');ax.legend(loc='lower left',handlelength=.8);data['panels']['b']={'filtered_endpoints':x,'excess':y}
    ax=fig.add_subplot(gs[1,1]);panel(ax,'c','How much GU repair recovers');x=[r['original_excess'] for r in cases];y=[r['fraction_excess_recovered'] for r in cases]
    for c,marker,col in zip(('incell','cellfree','kasugamycin'),('o','s','^'),(TEAL,BLUE,ORANGE)):
        rr=[r for r in cases if r['condition']==c];ax.scatter([r['original_excess'] for r in rr],[r['fraction_excess_recovered'] for r in rr],s=20,c=col,marker=marker)
    ax.set(xlabel=r'Original excess $E_{0.5}$',ylabel='Fraction recovered by GU addbacks',ylim=(-.04,1.04));data['panels']['c']={'original_excess':x,'GU_recovered_fraction':y};save(fig,'Figure_S3',out,data)

def figure6(out):
    result=load('results/expansion/barrier_pilot/validation.json');rows=result['rows'];summary=result['summary'];fig,axes=grid(height=158);data={'claim':'A prospective synthetic pilot fails engineering expansion gates and cannot distinguish saturated additive fields','scope':'24 new synthetic scaffolds;72 landscapes; no biological validation','panels':{},'summary':summary}
    ax=axes[0,0];panel(ax,'a','Frozen predictions versus folding')
    for variant,col,mark in (('baseline',GREY,'o'),('targeted',TEAL,'s'),('control',ORANGE,'^')):
        x=[r['predictions'][variant]['filtration']['longest_lifetime'] for r in rows];y=[r['validation'][variant]['filtration']['longest_lifetime'] for r in rows];ax.scatter(x,y,s=20,c=col,marker=mark,lw=.6,label=variant.capitalize());data['panels'][f'a_{variant}']={'predicted':x,'observed':y}
    ax.plot([0,.2],[0,.2],ls=':',c=GREY,lw=.7);ax.set_xscale('symlog',linthresh=.0001);ax.set_yscale('symlog',linthresh=.0001);ax.set(xlabel='Predicted lifetime (symlog)',ylabel='Ensemble lifetime\n(symlog)',xlim=(-.00002,.25),ylim=(-.00002,.25));ax.legend(loc='upper right',fontsize=7.5,handlelength=.8,labelspacing=.2,frameon=True,facecolor='white',edgecolor='white',framealpha=1.)
    ax=axes[0,1];panel(ax,'b','Targeted minus matched control');actual=[];predicted=[]
    for r in rows:
        actual.append(r['validation']['targeted']['filtration']['longest_lifetime']-r['validation']['control']['filtration']['longest_lifetime']);predicted.append(r['predictions']['targeted']['filtration']['longest_lifetime']-r['predictions']['control']['filtration']['longest_lifetime'])
    ax.plot(range(1,25),actual,'o',c=TEAL,ms=3,label='Ensemble');ax.plot(range(1,25),predicted,'x',c=ORANGE,ms=3,label='Prediction');ax.axhline(0,c=GREY,ls=':',lw=.7);ax.set_yscale('symlog',linthresh=.0001);ax.set(xlabel='Frozen scaffold order',ylabel=r'$L_{\rm target}-L_{\rm control}$'+'\n(symlog)',xticks=[1,8,16,24]);ax.legend(loc='lower left',bbox_to_anchor=(0,1.005),ncol=2,handlelength=.5,fontsize=7,borderaxespad=0,handletextpad=.3,columnspacing=.6);data['panels']['b']={'scaffolds':[r['scaffold'] for r in rows],'ensemble':actual,'prediction':predicted,'interval':summary['targeted_minus_control_lifetime']}
    axes[0,1].get_legend().remove()  # Encodings are defined in the manuscript caption.
    ax=axes[1,0];panel(ax,'c','Intervention direction and coverage');counts=[]
    for v in ('targeted','control'):
        r=[a for q in rows for a in q['intervention_outcomes'] if a['variant']==v];success=sum(a['covered'] and a['success'] for a in r);failed=sum(a['covered'] and not a['success'] for a in r);uncovered=sum(not a['covered'] for a in r);counts.append([success,failed,uncovered])
    matrix=np.array(counts);left=np.zeros(2)
    for k,col,label in zip(range(3),(TEAL,ORANGE,LIGHT),('Correct','Failed','Uncovered')):ax.barh([0,1],matrix[:,k],left=left,color=col,height=.55,label=label);left+=matrix[:,k]
    ax.set(yticks=[0,1],yticklabels=['Targeted','Control'],xlabel='Interventions (24 per class)',xlim=(0,24),ylim=(-.6,1.6));ax.invert_yaxis();ax.legend(loc='lower left',bbox_to_anchor=(-.55,-.7),ncol=3,fontsize=7.5,handlelength=.8,columnspacing=.7);data['panels']['c']={'counts_correct_failed_uncovered':counts}
    ax=axes[1,1];panel(ax,'d','Lifetime prediction error');x=[];z=[]
    for r in rows:
        x.append(np.mean([abs(r['predictions'][v]['filtration']['longest_lifetime']-r['validation'][v]['filtration']['longest_lifetime']) for v in ('baseline','targeted','control')]));z.append(np.mean([r['validation'][v]['filtration']['longest_lifetime'] for v in ('baseline','targeted','control')]))
    for i,(a,b) in enumerate(zip(x,z)):ax.plot([0,1],[a,b],c=LIGHT,lw=.7);ax.plot([0,1],[a,b],'o',c=TEAL,ms=2)
    ax.plot([0,1],[np.mean(x),np.mean(z)],'D',c=INK,ms=5);ax.set_yscale('symlog',linthresh=.0001);ax.set(xticks=[0,1],xticklabels=['Hairpin','Zero barrier'],ylabel='Scaffold mean absolute\nerror (symlog)',xlim=(-.3,1.3));data['panels']['d']={'hairpin_MAE_by_scaffold':list(map(float,x)),'zero_MAE_by_scaffold':list(map(float,z))}
    ax=axes[2,0];panel(ax,'e','Barriers at lifetime ≥ 0.01');conf=summary['confusion'];matrix=np.array([[conf['TN'],conf['FP']],[conf['FN'],conf['TP']]]);ax.imshow(matrix,cmap=LinearSegmentedColormap.from_list('confusion',['#F3F5F5',TEAL]),vmin=0,vmax=72,aspect='auto')
    for i,j in itertools.product(range(2),repeat=2):ax.text(j,i,str(matrix[i,j]),ha='center',va='center',color='white' if matrix[i,j]>35 else INK,fontsize=11)
    ax.set(xticks=[0,1],xticklabels=['No','Yes'],yticks=[0,1],yticklabels=['No','Yes'],xlabel='Predicted barrier',ylabel='Ensemble barrier');ax.tick_params(length=0);data['panels']['e']={'counts':matrix.tolist(),'definition':'rows observed no/yes;columns predicted no/yes; all72 landscapes'}
    ax=axes[2,1];panel(ax,'f','Correlation-aware reference');x=[];y=[]
    for r in rows:
        for v in ('baseline','targeted','control'):f=r['validation'][v];x.append(f['filtration']['integrated_excess']);y.append(f['references']['roughness5percent_mean_I_E'])
    ax.scatter(x,y,s=13,c=TEAL);limit=max(x+y)*1.1;ax.plot([0,limit],[0,limit],ls=':',c=GREY,lw=.8);ax.set(xlabel=r'Ensemble $I_E$',ylabel='Exact roughness-reference\n'+r'mean $I_E$',xlim=(-.001,limit),ylim=(-.001,limit));ax.text(.04,.9,'72/72 additive fits saturated',transform=ax.transAxes,fontsize=7.5);data['panels']['f']={'ensemble':x,'roughness_reference':y}
    save(fig,'Figure_6',out,data)

def supplement4(out):
    stress=load('results/expansion/residual_stress.json');fig,axes=grid(2,2,135);data={'claim':'Focused sensitivity checks do not turn the heat-reference residual into evidence for a folding-specific mechanism','panels':{}}
    for c,k,marker in zip(CONDITIONS,(-.09,.09),('o','s')):
        s=stress['conditions'][c]
        for j,t in enumerate(('0.01','0.05','0.1')):interval(axes[0,0],s['heat_roughness_sensitivities'][t],j+k,COLORS[c],marker)
        for j,n in enumerate((1,5)):interval(axes[0,1],s['heat_influence'][f'exclude_top{n}_absolute_cluster_contributions'],j+k,COLORS[c],marker)
        for j,kind in enumerate(('exclude_saturated','residual_df_at_least10percent_n')):interval(axes[1,0],s['codon'][kind],j+k,COLORS[c],marker)
        for j,t in enumerate(('1e-12','1e-10','1e-08')):interval(axes[1,1],s['tie_scale']['contrasts'][t],j+k,COLORS[c],marker)
        data['panels'][c]=s
    titles=('Tighter roughness matching','Largest cluster contributions removed','Saturation exclusions','Codon-fit numerical tie scale')
    for ax,letter,title in zip(axes.flat,'abcd',titles):panel(ax,letter,title);ax.axvline(0,c=GREY,ls=':',lw=.8);ax.invert_yaxis();ax.xaxis.set_major_locator(MaxNLocator(3))
    axes[0,0].set(yticks=range(3),yticklabels=['≤1%','≤5%','≤10%'],xlabel=r'Heat-reference − observed $I_E$')
    axes[0,1].set(yticks=[0,1],yticklabels=['Top 1','Top 5'],xlabel=r'Heat-reference − observed $I_E$')
    axes[1,0].set(yticks=[0,1],yticklabels=['df > 0','df/n ≥ 0.1'],xlabel=r'Additive − observed $I_E$')
    axes[1,1].set(yticks=range(3),yticklabels=[r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$'],xlabel=r'Additive − observed $I_E$')
    axes[0,0].legend([Line2D([],[],marker='o',c=TEAL,ls=''),Line2D([],[],marker='s',c=BLUE,ls='')],['In cell','Cell free'],loc='lower left',bbox_to_anchor=(-.25,-.85),ncol=2,handlelength=.8)
    save(fig,'Figure_S4',out,data)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--figures',default='1,2,3,4,5');a=p.parse_args();out=Path(a.out).resolve()
    if out==ROOT or ROOT in out.parents:raise ValueError('Figure exports must remain outside computational archive')
    out.mkdir(parents=True,exist_ok=True);(ROOT/'results/figure_data').mkdir(exist_ok=True)
    methods={'1':figure1,'2':figure2,'3':figure3,'4':figure4,'5':figure5,'6':figure6,'S1':supplement1,'S2':supplement2,'S3':supplement3,'S4':supplement4}
    for number in a.figures.split(','):SEEN.clear();methods[number](out)
