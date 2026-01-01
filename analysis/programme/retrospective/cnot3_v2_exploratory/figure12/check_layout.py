"""Inspect final Figure12 text geometry and deterministic export replay."""
import sys
if not __debug__:
    raise RuntimeError('Optimization refused')
sys.dont_write_bytecode=True
import hashlib
import importlib.util
import json
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib

HERE=Path(__file__).resolve().parent
DATA='3f40c1d627effd1fb0934b6ff676952bcca7030ae88e00b774b5576d061dbf31'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    qa=HERE/'qa_v5'
    if qa.exists():raise FileExistsError(qa)
    qa.mkdir()
    spec=importlib.util.spec_from_file_location('figure12_reviewed_renderer',HERE/'render_figure12.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    captures=[]
    old_close=plt.close
    plt.close=lambda value=None:captures.append(value) if hasattr(value,'findobj') else None
    sys.argv=[str(HERE/'render_figure12.py'),'--data',str(HERE/'panel_data/Figure_12_v5.json'),
              '--expected-data-sha256',DATA,'--out',str(qa/'replay')]
    module.main()
    plt.close=old_close
    fig=captures[-1];fig.canvas.draw(); renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    collisions=[]
    for i,t in enumerate(texts):
        a=t.get_window_extent(renderer)
        for u in texts[i+1:]:
            b=u.get_window_extent(renderer)
            dx=min(a.x1,b.x1)-max(a.x0,b.x0)
            dy=min(a.y1,b.y1)-max(a.y0,b.y0)
            if dx>.6 and dy>.6:
                collisions.append({'text_1':t.get_text(),'text_2':u.get_text(),'overlap_px':[dx,dy]})
    equal={name:sha(HERE/'final_v5'/name)==sha(qa/'replay'/name) for name in ['Figure_12.pdf','Figure_12.svg','Figure_12.png']}
    image=Image.open(HERE/'final_v5/Figure_12.png')
    image.convert('L').resize((1234,1163)).save(qa/'Figure_12_grayscale.png')
    receipt={'status':'passed' if not collisions and all(equal.values()) else 'failed',
       'text_pair_bbox_overlaps_greater_than_0_6_px':collisions,'visible_text_count':len(texts),
       'deterministic_replay_asset_equal':equal,'panel_data_sha256':sha(HERE/'panel_data/Figure_12_v5.json'),
       'PNG_dimensions_px':list(image.size),'PDF_physical_mm':[174,164],
       'minimum_visible_text_pt':min(t.get_fontsize() for t in texts),
       'panel_data_no_scientific_changes_from_layout_v1':None}
    original=json.loads((HERE/'Figure_12.json').read_text(encoding='utf-8'))
    final=json.loads((HERE/'panel_data/Figure_12_v5.json').read_text(encoding='utf-8'))
    receipt['panel_data_no_scientific_changes_from_layout_v1']=original['panels']==final['panels']
    if not receipt['panel_data_no_scientific_changes_from_layout_v1']:raise ValueError('Scientific panel data changed')
    (qa/'layout_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    old_close(fig)
    print(json.dumps(receipt,indent=2))
    if receipt['status']!='passed':raise RuntimeError('Layout/export check failed')
if __name__=='__main__':main()
