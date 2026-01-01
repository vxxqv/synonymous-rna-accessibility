from pathlib import Path
import json,csv
ROOT=Path(__file__).resolve().parents[1]
blocks=json.loads((ROOT/'data/native_windows.json').read_text())
transcripts={r['transcript_id']:r for r in json.loads((ROOT/'data/provenance/transcript_mapping.json').read_text())}
cds={(r['transcript_id'],r['locus_tag']):r for r in csv.DictReader((ROOT/'data/provenance/cds_mapping.csv').open())}
rows=[]
for b in blocks:
    tr=transcripts[b['transcript_id']];source=cds[b['transcript_id'],b['locus_tag']]
    partner={}
    for c,pairs in b['touching_pairs'].items():
        partner[c]={i:j for i,j in pairs};partner[c].update({j:i for i,j in pairs})
    for k,base in enumerate(b['native']):
        position=b['transcript_pos0']+k;context_position=b['variable_start0']+k
        row={'species':'Escherichia coli','strain':'K-12 MG1655','assembly':'U00096.2','annotation_date':'2013-02-26',
             'transcript_id':b['transcript_id'],'gene':b['gene'],'locus_tag':b['locus_tag'],'protein_id':source['protein_id'],
             'chromosome':'U00096.2','genomic_position1':tr['start']+position if b['strand']=='+' else tr['end']-position,
             'strand':b['strand'],'transcript_position1':position+1,'codon_index1':b['codon_index1']+k//3,'codon_phase1':k%3+1,
             'block_id':b['id'],'base':base,'mapping_confidence':'unique_exact_full_transcript','structural_evidence':'mirror-derived noncrossing target; original CT identity pending'}
        for c in ('incell','cellfree','kasugamycin'):
            p=partner[c].get(context_position)
            row[c+'_partner_transcript_position1']=None if p is None else b['context_start0']+p+1
            row[c+'_reactivity']=b['reactivity'][c][context_position]
        rows.append(row)
with (ROOT/'data/native_position_map.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print('mapped native variable positions',len(rows))
