"""Figure 11 from the immutable PERSIST prediction artifacts only.

No folding, source-table access, model fitting, exclusions or uncertainty
estimation occurs here. --repo-root accepts a reconstructed repository layout.
All graphical exports must be outside that computational repository.
"""
from pathlib import Path
import argparse
import hashlib
from itertools import combinations
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

SOURCES = {
    'analysis/programme/retrospective/persist_v2/protocol.json':
        'a3dc686a23acaa12af2f2f44154b8baac819cb6df3a0a6d81328dd3ebdc5fd24',
    'results/programme/persist_v2/train/features.json':
        '3d2f6cbc2804e4a07932f9850e6955c5dd2b0513dfcaa7a5f74f47ba0009229f',
    'results/programme/persist_v2/train/predictions.json':
        '2a22b0977afc031bfd2eceba15f5769daffad51f01a11319e81f787c7090a103',
    'results/programme/persist_v2/validate/validation.json':
        'a8aff2366b9735002cc6fb1c0ca844f3b04fcb009c2563acb79a4af93d98a86e'}
MODELS = ('development_median', 'development_mean', 'composition_ridge',
          'codon_composition_ridge', 'nearest_CDS', 'codon_composition_MUP_ridge')
LABELS = ('Median', 'Mean', 'Composition', 'Codon + comp.', 'Nearest CDS', 'Codon + MUP')
PRIMARY, COMPOSITION = MODELS[-1], MODELS[2]
INK, TEAL, BLUE, GREY, LIGHT, ORANGE = '#20272B', '#176B68', '#42699B', '#737C84', '#DCE1E4', '#B85E31'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(a, b, label):
    if not np.allclose(a, b, rtol=2e-13, atol=2e-15):
        raise ValueError('Saved numerical mismatch: '+label)


def read_frozen(repo):
    records = {}
    for relative, expected in SOURCES.items():
        if sha(repo/relative) != expected:
            raise ValueError('Frozen numerical source changed: '+relative)
        records[relative] = json.loads((repo/relative).read_text(encoding='utf-8'))
    protocol, features, training, validation = [records[key] for key in SOURCES]
    protocol_sha = SOURCES[next(iter(SOURCES))]
    if any(record['protocol_sha256'] != protocol_sha for record in (features,training,validation)):
        raise ValueError('Protocol identity differs between saved sources')
    if training['features_sha256'] != SOURCES['results/programme/persist_v2/train/features.json']:
        raise ValueError('Training feature identity differs')
    if validation['prediction_freeze']['prediction_sha256'] != SOURCES['results/programme/persist_v2/train/predictions.json']:
        raise ValueError('Validation prediction freeze differs')
    dev, held = protocol['development_ids'], protocol['holdout_ids']
    if len(dev) != 22 or len(held) != 15 or len(set(dev+held)) != 37:
        raise ValueError('Frozen 22/15 membership differs')
    if protocol['all_ids'] != dev+held or features['ids'] != dev+held or held != sorted(held):
        raise ValueError('Fixed identifier ordering differs')
    if training['development_ids'] != dev or training['holdout_ids'] != held or training['models'] != list(MODELS):
        raise ValueError('Frozen model or training membership differs')
    if [r['SequenceID'] for r in validation['rows']] != held:
        raise ValueError('Validation row ordering differs')
    if features['names'] != protocol['feature_names'] or len(features['names']) != 83:
        raise ValueError('Frozen feature columns differ')
    if protocol['primary_model'] != PRIMARY or protocol['models'] != list(MODELS):
        raise ValueError('Frozen comparison differs')
    sequence_hashes = {r['SequenceID']:r['normalized_full_sha256'] for r in protocol['cohort']}
    mup = {}
    for key in dev+held:
        vector = np.asarray(features['values'][key],dtype=float)
        report = features['diagnostics'][key]
        if vector.shape != (83,) or not np.all(np.isfinite(vector)):
            raise ValueError('Required feature row is invalid')
        value = float(vector[-1])
        close(value,report['value'],key+' MUP')
        close(value,report['row_formula_value'],key+' MUP alternative formula')
        if not 0 <= value <= 1 or report['normalized_sequence_sha256'] != sequence_hashes[key] or report['length'] != 928 or report['clipped']:
            raise ValueError('MUP input identity or numerical contract differs')
        mup[key] = value
    observed = np.asarray([r['target'] for r in validation['rows']],dtype=float)
    if not np.all(np.isfinite(observed)):
        raise ValueError('Required held-out values are nonfinite')
    maes, errors, predictions = {}, {}, {}
    for model in MODELS:
        prediction = np.asarray(training['predictions'][model],dtype=float)
        if prediction.shape != (15,) or not np.all(np.isfinite(prediction)):
            raise ValueError('Incomplete fixed prediction vector')
        close(prediction,[r['predictions'][model] for r in validation['rows']],model+' predictions')
        error = np.abs(prediction-observed)
        close(error,[r['absolute_errors'][model] for r in validation['rows']],model+' errors')
        mae = math.fsum(map(float,error))/15
        rmse = math.sqrt(math.fsum(float(v*v) for v in prediction-observed)/15)
        close(mae,validation['models'][model]['mae'],model+' MAE')
        close(rmse,validation['models'][model]['rmse'],model+' RMSE')
        maes[model],errors[model],predictions[model] = mae,error.tolist(),prediction.tolist()
    best = min(MODELS[:-1],key=lambda model:maes[model])
    hurdle = .9*maes[best]
    passes = all(maes[m] > 0 and maes[PRIMARY] <= .9*maes[m] for m in MODELS[:-1])
    if passes != validation['at_least_ten_percent_better_than_every_baseline'] or passes:
        raise ValueError('Failure framing differs from the frozen decision')
    if best != COMPOSITION:
        raise ValueError('Prespecified display of strongest comparator no longer matches source')
    for model in MODELS[:-1]:
        close(1-maes[PRIMARY]/maes[model],validation['fractional_gains'][model],model+' relative gain')
    difference = np.asarray(errors[PRIMARY])-np.asarray(errors[COMPOSITION])
    close(difference.mean(),maes[PRIMARY]-maes[COMPOSITION],'mean paired error difference')
    return {'panels':{
        'a':{'development_ids':dev,'holdout_ids':held,'MUP_by_id':mup,
             'horizontal_position':'fixed evenly spaced offsets permuted by identifier index (stride 9 for development, 7 for holdout); position within source has no scientific scale'},
        'b':{'model_order':list(MODELS),'display_labels':dict(zip(MODELS,LABELS)),'MAE':maes,
             'best_baseline':best,'ten_percent_hurdle':hurdle,'passes':passes},
        'c':{'holdout_ids':held,'observed':observed.tolist(),
             'primary_predictions':predictions[PRIMARY],'composition_predictions':predictions[COMPOSITION],
             'identity_line_only':True,'fitted_line':False},
        'd':{'holdout_ids':held,'primary_absolute_errors':errors[PRIMARY],
             'composition_absolute_errors':errors[COMPOSITION],
             'primary_minus_composition_absolute_error':difference.tolist(),
             'positive_means':'primary has larger absolute error'}},
        'additional_verified_summaries':{'primary_excess_over_composition':maes[PRIMARY]/maes[COMPOSITION]-1,
             'MUP_gain_over_matched_codon_composition':1-maes[PRIMARY]/maes[MODELS[3]]},
        'scope':'One selected reporter/context within one study. Source transfer, not independent genes. Absolute coefficient calibration unresolved.',
        'measurement':'Deposited cellular coefficient in original numeric units; MUP is a modeled ensemble expectation, not measured structure.',
        'verification':{'all37_MUP_values_checked':True,'all90_predictions_and_errors_checked':True,
             'all6_MAEs_and_RMSEs_recomputed':True,'frozen_gain_decision_checked':True,
             'original_source_table_opened':False,'folding':False,'model_refitting':False,
             'exclusions':False,'inferential_intervals':False}}


def inspect_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    clipping, overlaps, small = [], [], []
    for text in texts:
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.x1+1 or box.y1 > fig.bbox.y1+1:
            clipping.append(text.get_text())
        if text.get_fontsize() < 8:
            small.append(text.get_text())
    for a,b in combinations(texts,2):
        x,y = a.get_window_extent(renderer),b.get_window_extent(renderer)
        if min(x.x1,y.x1)-max(x.x0,y.x0) > .5 and min(x.y1,y.y1)-max(x.y0,y.y0) > .5:
            overlaps.append([a.get_text(),b.get_text()])
    partial_alpha = []
    for item in fig.findobj():
        if hasattr(item,'get_alpha') and item.get_alpha() not in (None,0,1):
            partial_alpha.append(type(item).__name__)
    return {'clipped_text':clipping,'overlapping_text':overlaps,'font_below_8pt':small,
            'partial_alpha_artists':partial_alpha,'visible_text_items':len(texts)}


def render(repo,out,qa,manifest):
    repo,out,qa,manifest = [Path(p).resolve() for p in (repo,out,qa,manifest)]
    if out.is_relative_to(repo) or qa.is_relative_to(repo):
        raise ValueError('Graphical exports and QA must be outside the computational repository')
    for path in (out,qa,manifest.parent):
        path.mkdir(parents=True,exist_ok=True)
    data = read_frozen(repo)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
        'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,
        'axes.linewidth':.65,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'PERSIST-v2-Figure-11','savefig.facecolor':'white',
        'legend.fontsize':8,'legend.frameon':False,'legend.framealpha':1,'axes.unicode_minus':False})
    fig = plt.figure(figsize=(174/25.4,158/25.4),facecolor='white')
    a = fig.add_axes([.125,.645,.34,.265])
    b = fig.add_axes([.70,.645,.27,.265])
    c = fig.add_axes([.12,.135,.29,.32])
    d = fig.add_axes([.70,.135,.27,.32])
    for x,y,letter,title in ((.035,.957,'a','Modeled structural input'),(.55,.957,'b','Frozen prediction test'),
                             (.035,.540,'c','Held-out predictions'),(.55,.540,'d','Paired error differences')):
        fig.text(x,y,letter,fontsize=10,weight='bold',va='center')
        fig.text(x+.03,y,title,fontsize=8.5,weight='bold',va='center')

    pa = data['panels']['a']
    for center,ids,color,marker in ((0,pa['development_ids'],GREY,'o'),(1,pa['holdout_ids'],BLUE,'D')):
        stride = 9 if len(ids)==22 else 7
        offsets = np.linspace(-.24,.24,len(ids))[(np.arange(len(ids))*stride)%len(ids)]
        pa.setdefault('display_x_by_id',{}).update({key:float(center+x) for key,x in zip(ids,offsets)})
        a.scatter(center+offsets,[pa['MUP_by_id'][key] for key in ids],s=22,marker=marker,
                  facecolor=color,edgecolor='white',linewidth=.35,zorder=3)
    a.set(xlim=(-.45,1.45),ylim=(.24,.47),yticks=[.25,.30,.35,.40,.45],
          xticks=[0,1],xticklabels=['Development\n22 constructs','Held-out\n15 constructs'],
          ylabel='Mean unpaired probability')
    a.tick_params(axis='x',length=0,pad=6)

    pb = data['panels']['b']
    vals = [pb['MAE'][name] for name in MODELS]
    colors = [TEAL if name==PRIMARY else BLUE if name==COMPOSITION else GREY for name in MODELS]
    b.barh(range(6),vals,height=.60,color=colors,zorder=2)
    b.set_yticks(range(6),LABELS)
    b.invert_yaxis()
    b.set(xlim=(0,.058),xticks=[0,.02,.04],xlabel='Held-out MAE')
    for k,(name,value) in enumerate(zip(MODELS,vals)):
        b.text(value+.0012,k,f'{value:.5f}',va='center',fontsize=8,
               color=colors[k],weight='bold' if name==PRIMARY else 'normal')
    b.axvline(pb['ten_percent_hurdle'],color=ORANGE,lw=1,ls=(0,(3,2)),zorder=3)
    fig.text(.70,.925,f'10% hurdle: {pb["ten_percent_hurdle"]:.5f}',color=ORANGE,fontsize=8)
    for label in b.get_yticklabels():
        if label.get_text()==LABELS[-1]: label.set_color(TEAL); label.set_weight('bold')

    pc = data['panels']['c']
    c.plot([.08,.215],[.08,.215],color=GREY,lw=.9,ls=(0,(3,2)),zorder=1)
    for x,y1,y2 in zip(pc['observed'],pc['primary_predictions'],pc['composition_predictions']):
        c.plot([x,x],[y1,y2],color=LIGHT,lw=.8,zorder=2)
    c.scatter(pc['observed'],pc['composition_predictions'],s=32,marker='D',facecolors='none',
              edgecolors=BLUE,linewidths=.9,zorder=4)
    c.scatter(pc['observed'],pc['primary_predictions'],s=20,marker='o',facecolors=TEAL,
              edgecolors='white',linewidths=.25,zorder=5)
    c.set(xlim=(.08,.215),ylim=(.08,.215),xticks=[.10,.15,.20],yticks=[.10,.15,.20],
          xlabel='Observed coefficient',ylabel='Predicted coefficient',aspect='equal')
    c.legend(handles=[Line2D([],[],marker='o',ls='none',ms=4.5,mfc=TEAL,mec='white',label='Primary'),
                      Line2D([],[],marker='D',ls='none',ms=4.5,mfc='none',mec=BLUE,label='Composition'),
                      Line2D([],[],ls=(0,(3,2)),color=GREY,lw=.9,label='Identity')],
             loc='upper left',borderaxespad=.25,handlelength=1.3,handletextpad=.4,labelspacing=.22)

    pd = data['panels']['d']
    diff = np.asarray(pd['primary_minus_composition_absolute_error'])
    color = [ORANGE if v>0 else TEAL if v<0 else GREY for v in diff]
    d.barh(range(15),diff,color=color,height=.58,zorder=2)
    d.axvline(0,color=INK,lw=.8,zorder=3)
    d.set_yticks(range(15),pd['holdout_ids'])
    d.invert_yaxis()
    d.set(xlim=(-.008,.0102),xticks=[-.005,0,.005,.010],
          xticklabels=['-0.005','0','0.005','0.010'],
          xlabel='Absolute error difference\n(primary minus composition)')
    d.tick_params(axis='y',length=0,pad=4)
    fig.text(.70,.491,'Lower error',color=TEAL,fontsize=8)
    fig.text(.85,.491,'Higher error',color=ORANGE,fontsize=8)
    fig.text(.035,.023,'One selected reporter cohort. Coefficient calibration unresolved; MUP is modeled, not measured structure.',fontsize=8)

    qa_report = inspect_layout(fig)
    fig.savefig(qa/'Figure_11_preview.png',dpi=180)
    (qa/'layout_check.json').write_text(json.dumps(qa_report,indent=2)+'\n',encoding='utf-8')
    if any(qa_report[key] for key in ('clipped_text','overlapping_text','font_below_8pt','partial_alpha_artists')):
        raise ValueError('Layout or alpha check failed: '+json.dumps(qa_report))
    for extension in ('pdf','svg','eps','png'):
        metadata = {'CreationDate':None,'ModDate':None} if extension=='pdf' else {'Date':None} if extension=='svg' else None
        fig.savefig(out/('Figure_11.'+extension),dpi=610,metadata=metadata,facecolor='white',transparent=False)
        if extension=='eps':
            # Remove the automatically added wall-clock comment only. Graphics
            # are unchanged; repeated exports do not encode render time.
            path = out/'Figure_11.eps'
            path.write_bytes(b''.join(line for line in path.read_bytes().splitlines(keepends=True)
                                      if not line.startswith(b'%%CreationDate:')))
    for relative,expected in SOURCES.items():
        if sha(repo/relative) != expected:
            raise RuntimeError('Numerical source changed during rendering')
    data['provenance'] = {'source_sha256':SOURCES,'renderer_sha256':sha(__file__),
                          'matplotlib':matplotlib.__version__,'numpy':np.__version__}
    data['export'] = {'width_mm':174,'height_mm':158,'minimum_ordinary_font_pt':8,'png_dpi':610,
                      'caption_external':True,'date_metadata_removed':True,'stable_svg_ids':True,'layout_check':qa_report,
                      'sha256':{ext:sha(out/('Figure_11.'+ext)) for ext in ('pdf','svg','eps','png')}}
    manifest.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (qa/'render_receipt.json').write_text(json.dumps({'manifest_sha256':sha(manifest),
        'provenance':data['provenance'],'export':data['export'],'verification':data['verification']},indent=2)+'\n',encoding='utf-8')
    plt.close(fig)
    print(json.dumps({'figure':str(out/'Figure_11.pdf'),'manifest':str(manifest),
                      'layout':qa_report,'comparison':pb,'real_analysis_performed':False}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--qa',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    arguments = parser.parse_args()
    render(arguments.repo_root,arguments.outdir,arguments.qa,arguments.manifest)
