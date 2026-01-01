"""Read only explicitly whitelisted metadata/sequence cells, never phenotype cells."""
from pathlib import Path
import json,csv,collections,openpyxl,hashlib

if not __debug__:raise RuntimeError('Optimized Python prohibited')
p=Path('work/research_cycle8/cprime_gate')

def selected_xlsx(file,sheet,columns,start=2):
    w=openpyxl.load_workbook(p/file,read_only=True,data_only=False)
    s=w[sheet]
    # Read each approved column separately. No phenotype-column cell is returned.
    vals={name:[c[0] for c in s.iter_rows(min_row=start,min_col=col,max_col=col,values_only=True)] for name,col in columns.items()}
    n=max(map(len,vals.values()),default=0)
    rows=[{k:(v[i] if i<len(v) else None) for k,v in vals.items()} for i in range(n)]
    return [r for r in rows if any(x is not None for x in r.values())]

allowed={'combined_var_id':1,'var_id':2,'hgvsc':3,'hgvsp':4,'exon_number':5,'syn_var':6,'cell_line':7,'consequence':27,'Target_sequence':33}
multi=selected_xlsx('supp_Supplementary Data 8.xlsx','Supplementary_Table_3_Multi_exo',allowed)
notation={}
for n in (7,8):
    notation[str(n)]=selected_xlsx('supp_Supplementary Data '+str(n)+'.xlsx','cDNA notation',{'combined_var_id':1,'transcript':2,'manuscript_id':3,'variant_id_edit':4,'variant_id_marker':5,'hgvsp':6,'hgvsc_edit':7,'hgvsc_marker':8})

validation={}
for sheet,col in [('technical valid KBM7 K562 Fc',3),('validation KBM7',4),('validation K562',3)]:
    validation[sheet]=selected_xlsx('supp_Supplementary Data 6.xlsx',sheet,{'Paper_name':col})
# The names sheet contains metadata only at these explicit columns; exclude codon numeric cells.
name_map=selected_xlsx('supp_Supplementary Data 6.xlsx','names',{'source_id':1,'paper_label':3,'display_id':5},start=1)

contexts=collections.defaultdict(set)
for r in multi:
    if r['combined_var_id'] and r['cell_line']:contexts[str(r['cell_line'])].add(r['combined_var_id'])
design=json.loads((p/'multi_exon_selected_metadata.json').read_text())['rows']
groups=collections.defaultdict(set)
for r in design:groups[r['var_id']].add(r['combined_var_id'])
summary={'design_rows':len(design),'intended_edits':len(groups),'unique_combined_genotypes':len({r['combined_var_id'] for r in design}),
         'intended_edit_combination_multiplicity':dict(collections.Counter(len(v) for v in groups.values())),
         'multi_source_metadata_rows':len(multi),'context_unique_combined_ids':{k:len(v) for k,v in contexts.items()},
         'two_context_combined_id_intersection':len(set.intersection(*contexts.values())) if len(contexts)==2 else None,
         'context_names':sorted(contexts),'annotation_rows':{k:len(v) for k,v in notation.items()},
         'no_phenotype_values_or_missingness_inspected':True,'no_complete_haplotype_or_isoform_claim':True}
payload={'selected_columns':{'multi_data8':allowed,'data6_ID_only':validation.keys().__iter__().__class__.__name__},
         'multi_metadata':multi,'notation_metadata':notation,'validation_ID_columns':validation,'name_map':name_map,
         'summary':summary}
payload['selected_columns']['data6_ID_only']={'technical valid KBM7 K562 Fc':3,'validation KBM7':4,'validation K562':3,'names':[1,3,5]}
(p/'selected_assay_metadata.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
(p/'metadata_contract_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
print('VALIDATION_ID_ONLY',json.dumps(validation,indent=2))
print('NAME_MAP',json.dumps(name_map,indent=2))
