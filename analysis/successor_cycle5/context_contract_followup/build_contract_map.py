"""Build an assay-accession map from metadata only, with unresolved units null."""
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'mordstein_eligibility'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def attributes(path):
    root=ET.parse(path).getroot()
    return {n.findtext('TAG'):n.findtext('VALUE') for n in root.iter('SAMPLE_ATTRIBUTE')}

def main():
    runs=json.loads((HERE/'ena_runs.json').read_text())
    historical=HERE.parent/'splicing_sources/mordstein_sra_metadata.response'
    old=list(csv.DictReader(historical.read_text().splitlines(),delimiter='\t'))
    common=('run_accession','sample_accession','library_name','experiment_title')
    if {tuple(r[k] for k in common) for r in runs}!={tuple(r[k] for k in common) for r in old}:
        raise RuntimeError('Current and prior run metadata differ')
    fractionation=sorted((r for r in runs if r['library_name'].startswith('CM6')),key=lambda r:r['library_name'])
    mapping=[]
    for row in fractionation:
        sample=row['sample_accession']
        source=HERE/(sample+'.xml')
        if not source.exists():source=PRIOR/(sample+'.xml')
        a=attributes(source)
        experiment=ET.parse(HERE/(row['experiment_accession']+'.xml')).getroot()
        lib=experiment.findtext('.//LIBRARY_NAME')
        if lib!=row['library_name']:raise RuntimeError('Library identity mismatch')
        sample_alias=experiment.findtext('.//SAMPLE_DESCRIPTOR/IDENTIFIERS/SUBMITTER_ID')
        mapping.append({**row,'structured_cell_line':a['cell_line'],'structured_assay':a['experiment type'],
            'intron':a['5UTR_intron'],'fraction':a['sample_type'],'sample_alias':sample_alias,
            'source_sample_file':str(source.relative_to(HERE.parent)),'source_sample_sha256':sha(source),
            'source_experiment_sha256':sha(HERE/(row['experiment_accession']+'.xml')),
            'label_suffix':lib[-1],'independent_preparation_ID':None,'matched_preparation_group_ID':None,
            'exact_variant_CDS_index':None,'prepared_same_culture_across_fractions':None,
            'title_vs_structured_cell_line_conflict':'HEK293' in row['experiment_title'] and a['cell_line']=='HeLa'})
    contract={'status':'NO_GO_EXACT_MATCHED_UNIT_CONTRACT_UNRESOLVED',
        'study':'PRJNA596086','secondary_study':'SRP237851','run_count':len(runs),
        'biosample_count':len({r['sample_accession'] for r in runs}),
        'current_run_metadata_agrees_with_prior':True,'fractionation_run_count':len(mapping),
        'fractionation_experiment_count':len({r['experiment_accession'] for r in mapping}),
        'fractionation_biosample_count':len({r['sample_accession'] for r in mapping}),
        'fractionation_mapping':mapping,'complete_same_CDS_unit_count':None,'unique_CDS_count':None,
        'independently_authenticated_preparation_count':None,'nominal_reported_variants':217,
        'source_category_sum':219,'category_overlap_or_exclusions_resolved':False,
        'analysis_accession_query_count':len(json.loads((HERE/'ena_analyses.json').read_text())),
        'exact_assay_variant_intersections':None,'outcome_analysis_authorized':False,
        'unknown_counts_are_not_zero':True}
    with (HERE/'contract_map.json').open('x',encoding='utf-8') as f:json.dump(contract,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({k:v for k,v in contract.items() if k!='fractionation_mapping'},indent=2))

if __name__=='__main__':main()
