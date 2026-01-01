"""Figure 9 from frozen saved pilot vectors only, with independent value checks.

This renderer never opens the source workbook, refits a model, or selects cases.
All 32 genotypes, seven models, and nine dependent library contrasts are shown.
"""
from pathlib import Path
import argparse
import hashlib
from itertools import combinations, product
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('results/programme/gill_reference_pilot')
EXPECTED = {
    'training/predictions.json': 'ff04894c283310fc9f47ae01f35191bc6a9be83493e39baaaf29d6da9bdd3660',
    'training/receipt.json': 'b2877db7763163aef44d8a43e13ab75680040608f94cd4d20bad5fbfcbe9ce61',
    'validation/validation.json': 'eedb95ab88b1e9860ba304da9921da3eee282afe1bd2d30088e412109aac6766',
    'validation/joint_scenarios.json': '78876369e146cef823a1445d18aa65c5d057ab89cdf1eff05e9436e3900b14dd',
    'validation/receipt.json': 'ed65af6f4460ffb5a839642237546b9c0938ccc885a366f25bb7b01a2e970a77'}
INK, TEAL, GREY, LIGHT = '#20272B', '#176B68', '#737C84', '#DCE1E4'
BLUE, ORANGE = '#42699B', '#B85E31'
MODEL_LABELS = {'constant': 'Constant', 'gc': 'GC', 'additive': 'Additive',
                'isotonic_additive': 'Isotonic additive', 'pairwise_half': 'Pairwise half',
                'heat': 'Heat', 'training_mean': 'Training mean'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_close(observed, expected, description):
    if not np.allclose(observed, expected, atol=2e-13, rtol=2e-13):
        raise ValueError('Saved result mismatch: '+description)


def order_energies(field):
    """Tensor conditional-mean contrasts, recomputed independently of the runner."""
    tensor = np.asarray(field).reshape(2, 4, 4)
    components = {(): np.full_like(tensor, tensor.mean())}
    result = []
    for order in (1, 2, 3):
        energy = 0.0
        for subset in combinations(range(3), order):
            complement = tuple(i for i in range(3) if i not in subset)
            marginal = tensor.mean(axis=complement, keepdims=True) if complement else tensor.copy()
            component = np.broadcast_to(marginal, tensor.shape).copy()
            for smaller, part in components.items():
                if len(smaller) < len(subset) and set(smaller).issubset(subset):
                    component -= part
            components[subset] = component
            energy += float(np.mean(component**2))
        result.append(energy)
    verify_close(sum(result), float(np.var(tensor)), 'ANOVA energy sum')
    return result


def maximum_barrier(field, graph):
    """Recover capacities through critical-threshold connectivity, not Floyd."""
    scores = np.asarray(field, dtype=float)
    capacity = np.full((32, 32), -np.inf)
    for threshold in sorted(set(scores)):
        unseen = {i for i in range(32) if scores[i] >= threshold}
        while unseen:
            seed = min(unseen)
            todo, component = [seed], {seed}
            unseen.remove(seed)
            while todo:
                node = todo.pop()
                for neighbor in graph[node]:
                    if neighbor in unseen:
                        unseen.remove(neighbor)
                        component.add(neighbor)
                        todo.append(neighbor)
            indices = sorted(component)
            capacity[np.ix_(indices, indices)] = threshold
    return float(np.max(np.minimum(scores[:, None], scores[None, :])-capacity))


def load_data():
    files = {}
    for relative, expected in EXPECTED.items():
        path = ROOT/BASE/relative
        if sha(path) != expected:
            raise ValueError('Frozen figure source changed: '+relative)
        files[relative] = json.loads(path.read_text(encoding='utf-8'))
    training = files['training/predictions.json']
    validation = files['validation/validation.json']
    scenarios = files['validation/joint_scenarios.json']
    for prefix, receipt in (('training', files['training/receipt.json']), ('validation', files['validation/receipt.json'])):
        for name in ('predictions.json',) if prefix == 'training' else ('validation.json', 'joint_scenarios.json'):
            if receipt['outputs'][name] != EXPECTED[prefix+'/'+name]:
                raise ValueError('Source receipt output mismatch')
    if validation['prediction_sha256'] != EXPECTED['training/predictions.json'] or training['identity'] != validation['identity']:
        raise ValueError('Training and validation freeze identity mismatch')
    keys = training['keys']
    if len(keys) != 32 or keys != validation['keys'] or len(set(keys)) != 32:
        raise ValueError('Figure requires all 32 matching genotypes')
    reference = training['reference_index']
    all_rmse, criteria = {}, {}
    for c in ('0.25', '0.5', '1.0'):
        outcome = np.asarray(validation['results'][c]['test_anchored'])
        estimates = training['predictions'][c]['models']
        all_rmse[c] = {}
        if list(estimates) != list(MODEL_LABELS) or outcome.shape != (32,) or outcome[reference] != 0:
            raise ValueError('Incomplete frozen model or anchored outcome array')
        for model, field in estimates.items():
            error = np.asarray(field)-outcome
            verify_close(error, validation['results'][c]['models'][model]['errors'], model+' error vector')
            rmse = math.sqrt(math.fsum(float(x*x) for x in error)/32)
            verify_close(rmse, validation['results'][c]['models'][model]['rmse32'], model+' RMSE')
            all_rmse[c][model] = rmse
        best = min(value for model, value in all_rmse[c].items() if model != 'pairwise_half')
        target = all_rmse[c]['pairwise_half']
        success = best > 0 and target <= .9*best
        if success != validation['summary']['all_pseudocount_decisions'][c]['passes_10_percent_hurdle']:
            raise ValueError('Primary decision mismatch')
        criteria[c] = {'best_baseline': best, 'target': target, 'hurdle': .9*best,
                       'passes': success, 'relative_excess': target/best-1}
    if any(row['passes'] for row in criteria.values()):
        raise ValueError('Figure failure framing no longer matches saved criteria')
    train_field = training['predictions']['0.5']['training_anchored']
    test_field = validation['results']['0.5']['test_anchored']
    energy = {'training': order_energies(train_field), 'heldout': order_energies(test_field)}
    verify_close(energy['heldout'], [validation['results']['0.5']['observed_geometry']['anova_mean_square_by_order'][str(i)] for i in (1, 2, 3)], 'held-out ANOVA')
    graph = [[j for j, y in enumerate(keys) if sum(a != b for a, b in zip(x, y)) == 1] for x in keys]
    if set(map(len, graph)) != {7}:
        raise ValueError('Saved genotype graph differs from frozen design')
    barriers = []
    for rna in ('B', 'C', 'D'):
        row = []
        for dna in ('E', 'F', 'G'):
            scenario = scenarios['0.5'][rna+'/'+dna]
            value = maximum_barrier(scenario['field'], graph)
            verify_close(value, scenario['geometry']['maximum_barrier'], rna+'/'+dna+' maximum barrier')
            row.append(value)
        barriers.append(row)
    return {'keys': keys, 'reference_index': reference, 'reference_key': training['reference_key'],
            'all_rmse': all_rmse, 'criteria': criteria, 'primary_pseudocount': .5,
            'prediction': training['predictions']['0.5']['models']['pairwise_half'], 'observed': test_field,
            'energy': energy, 'scenario_barriers': barriers,
            'scientific_scope': validation['summary']['scope'], 'source_identity': training['identity']}


def render(out, qa):
    if not __debug__:
        raise RuntimeError('Integrity checks require normal Python execution')
    data = load_data()
    out, qa = Path(out), Path(qa)
    out.mkdir(parents=True, exist_ok=True)
    qa.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .65,
        'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'savefig.facecolor': 'white', 'legend.fontsize': 8, 'legend.frameon': False})
    fig = plt.figure(figsize=(174/25.4, 143/25.4))
    axa = fig.add_axes([.185, .625, .32, .275])
    axb = fig.add_axes([.665, .625, .30, .275])
    axc = fig.add_axes([.185, .17, .32, .275])
    axd = fig.add_axes([.665, .17, .24, .275])
    for x, y, letter, title in ((.035, .953, 'a', 'Frozen prediction test'),
                               (.585, .953, 'b', 'Held-out reproduction'),
                               (.035, .497, 'c', 'Orthogonal score components'),
                               (.585, .497, 'd', 'Dependent library sensitivity')):
        fig.text(x, y, letter, fontsize=10, weight='bold', va='center')
        fig.text(x+.03, y, title, fontsize=8.5, weight='bold', va='center')

    values = data['all_rmse']['0.5']
    labels = list(MODEL_LABELS)
    colors = [TEAL if name == 'pairwise_half' else INK if name == 'constant' else GREY for name in labels]
    axa.barh(range(7), [values[name] for name in labels], color=colors, height=.62)
    axa.set_yticks(range(7), [MODEL_LABELS[name] for name in labels])
    axa.invert_yaxis()
    axa.set_xlim(0, .79)
    axa.set_xticks([0, .2, .4, .6])
    axa.set_xlabel('RMSE over 32 genotypes (log$_2$)')
    for i, name in enumerate(labels):
        axa.text(values[name]+.015, i, f'{values[name]:.3f}', va='center', fontsize=8,
                 weight='bold' if name == 'pairwise_half' else 'normal', color=colors[i])
    hurdle = data['criteria']['0.5']['hurdle']
    axa.axvline(hurdle, color=ORANGE, ls=(0, (3, 2)), lw=1)
    axa.text(.98, 1.04, f'10% hurdle: {hurdle:.3f}', color=ORANGE,
             ha='right', transform=axa.transAxes, fontsize=8)
    for tick in axa.get_yticklabels():
        if tick.get_text() == MODEL_LABELS['pairwise_half']:
            tick.set_weight('bold')
            tick.set_color(TEAL)

    observed = np.asarray(data['observed'])
    predicted = np.asarray(data['prediction'])
    ref = data['reference_index']
    mask = np.arange(32) != ref
    axb.plot([-1, 1.1], [-1, 1.1], color=GREY, lw=.9, ls=(0, (3, 2)), zorder=1)
    axb.scatter(observed[mask], predicted[mask], s=21, facecolor=TEAL, edgecolor='white', linewidth=.35, zorder=3)
    axb.scatter([0], [0], marker='D', s=38, facecolor='white', edgecolor=INK, linewidth=.9, zorder=4)
    axb.set(xlim=(-1, 1.1), ylim=(-1, 1.1), xticks=[-1, 0, 1], yticks=[-1, 0, 1],
            xlabel='Observed contrast (log$_2$)', ylabel='Predicted contrast (log$_2$)')
    axb.legend(handles=[Line2D([], [], marker='D', ls='none', ms=5, mfc='white', mec=INK, label='Reference anchor'),
                        Line2D([], [], ls=(0, (3, 2)), color=GREY, lw=.9, label='Identity')],
               loc='upper left', borderaxespad=.25, handlelength=1.4, handletextpad=.5, labelspacing=.25)

    index = np.arange(3)
    width = .33
    axc.bar(index-width/2, data['energy']['training'], width, color=BLUE, label='Training')
    axc.bar(index+width/2, data['energy']['heldout'], width, facecolor='white', edgecolor=TEAL,
            hatch='///', linewidth=.9, label='Held-out')
    axc.set_xticks(index, ['Main', 'Pairwise', 'Three-way'])
    axc.set_ylim(0, .115)
    axc.set_yticks([0, .05, .10], ['0', '0.05', '0.10'])
    axc.set_xlabel('Factorial order')
    axc.set_ylabel('Mean square (log$_2^2$)')
    training_labels = index-width/2-.10
    training_labels[0] += .10
    for positions, field in ((training_labels, data['energy']['training']), (index+width/2+.10, data['energy']['heldout'])):
        for x, y in zip(positions, field):
            axc.text(x, y+.0035, f'{y:.3f}', ha='center', va='bottom', fontsize=8)
    axc.legend(loc='upper right', ncol=1, borderaxespad=.1, handlelength=1.2, handletextpad=.5, labelspacing=.25)

    matrix = np.asarray(data['scenario_barriers'])
    cmap = LinearSegmentedColormap.from_list('gill_barriers', ['#F3F6F5', '#7CAAA5', TEAL])
    im = axd.imshow(matrix, vmin=0, vmax=.5, cmap=cmap, interpolation='nearest', aspect='auto')
    axd.set_xticks(range(3), ['DNA 8', 'DNA 9', 'DNA 10'])
    axd.set_yticks(range(3), ['RNA 1', 'RNA 2', 'RNA 3'])
    axd.tick_params(length=0, pad=5)
    axd.set_xlabel('DNA denominator')
    for spine in axd.spines.values():
        spine.set_visible(False)
    for row, col in product(range(3), repeat=2):
        axd.text(col, row, f'{matrix[row, col]:.3f}', ha='center', va='center',
                 color='white' if matrix[row, col] > .33 else INK, fontsize=8.5)
    axd.add_patch(Rectangle((1.5, 1.5), 1, 1, fill=False, edgecolor=INK, linewidth=1.2))
    color_ax = fig.add_axes([.925, .17, .014, .275])
    colorbar = fig.colorbar(im, cax=color_ax, ticks=[0, .25, .5])
    colorbar.ax.set_yticklabels(['0', '0.25', '0.50'])
    colorbar.outline.set_linewidth(.5)

    fig.text(.185, .075, 'Raw interactions do not identify a mechanism.', fontsize=8, color=GREY)
    fig.text(.665, .075, 'Maximum barrier (log$_2$)', fontsize=8)
    fig.text(.035, .026, 'Same 32 genotypes; pseudocount 0.5. Shared libraries are dependent scenarios, not confidence intervals.', fontsize=8)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    clipping, overlaps = [], []
    for text in texts:
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.x1+1 or box.y1 > fig.bbox.y1+1:
            clipping.append(text.get_text())
        if text.get_fontsize() < 8:
            raise ValueError('Font below 8 pt')
    for a, b in combinations(texts, 2):
        x, y = a.get_window_extent(renderer), b.get_window_extent(renderer)
        if min(x.x1, y.x1)-max(x.x0, y.x0) > .5 and min(x.y1, y.y1)-max(x.y0, y.y0) > .5:
            overlaps.append([a.get_text(), b.get_text()])
    if clipping or overlaps:
        raise ValueError(json.dumps({'clipping': clipping, 'overlaps': overlaps}))
    for ext in ('pdf', 'svg', 'eps', 'png'):
        metadata = {'CreationDate': None, 'ModDate': None} if ext == 'pdf' else {'Date': None} if ext == 'svg' else None
        fig.savefig(out/('Figure_9.'+ext), dpi=610, metadata=metadata)
    fig.savefig(qa/'Figure_9_review.png', dpi=180)
    for i, ax in enumerate((axa, axb, axc, axd)):
        extent = ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.04, 1.07)
        fig.savefig(qa/f'panel_{chr(97+i)}.png', bbox_inches=extent, dpi=250)
    source_hashes = {str(BASE/relative).replace('\\', '/'): value for relative, value in EXPECTED.items()}
    if any(sha(ROOT/BASE/path) != value for path, value in EXPECTED.items()):
        raise RuntimeError('Source changed during render')
    data['provenance'] = {'renderer_sha256': sha(__file__), 'inputs': source_hashes,
                          'matplotlib': matplotlib.__version__, 'numpy': np.__version__}
    data['export'] = {'width_mm': 174, 'height_mm': 143, 'minimum_ordinary_font_pt': 8,
                      'clipped_text': clipping, 'text_overlap': overlaps,
                      'artifacts': {ext: sha(out/('Figure_9.'+ext)) for ext in ('pdf', 'svg', 'eps', 'png')}}
    data['verification'] = {'all_21_rmse_recomputed': True, 'all_672_errors_checked': True,
                            'energies_recomputed_by_conditional_means': True,
                            'nine_barriers_recomputed_by_threshold_connectivity': True,
                            'source_workbook_opened': False, 'model_refitting': False}
    destination = ROOT/'results/figure_data/Figure_9.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    (qa/'receipt.json').write_text(json.dumps({'data': str(destination), 'data_sha256': sha(destination),
        'provenance': data['provenance'], 'export': data['export'], 'verification': data['verification']}, indent=2)+'\n', encoding='utf-8')
    plt.close(fig)
    print(json.dumps({'figure': str(out/'Figure_9.pdf'), 'clipping': clipping, 'overlaps': overlaps,
                      'primary': data['criteria']['0.5']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--qa', type=Path, required=True)
    args = parser.parse_args()
    render(args.outdir, args.qa)
