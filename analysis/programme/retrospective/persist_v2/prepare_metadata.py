"""Prepare metadata and installed model settings only. No features or outcomes."""
from pathlib import Path
import json
import platform
import RNA
import numpy as np
from persist_prediction import sha256, text_sha, normalize_rna, snapshot_md, groups_at_distance, write_new

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = HERE.parents[2]
PROTEIN = 'dd0e8ca4f9a41c1fdd610ac0e48d953ff55be9483e2067db5e713804931fc148'
CONTEXT = 'e682ac5099bf1b6f4beb22f1f9817c89022c0bd31b279975a872546eddf726e9'


def main():
    if (HERE/'metadata_settings.json').exists() or (HERE/'Turner2004.par').exists():
        raise FileExistsError('Do not replace the model/metadata snapshot')
    eligibility = json.loads((BASE/'eligibility_results.json').read_text())
    families = [x for x in eligibility['families'] if x['protein_sha256'] == PROTEIN and x['context_sha256'] == CONTEXT]
    assert len(families) == 1
    ids = sorted(families[0]['longamp_ids'])
    development, holdout = [i for i in ids if i.startswith('310')], [i for i in ids if i.startswith('930')]
    assert len(ids) == 37 and len(development) == 22 and len(holdout) == 15 and sorted(development+holdout) == ids
    metadata = json.loads((BASE/'permitted_metadata_verified.json').read_text())
    lookup = {r['SequenceID']:r for r in metadata}
    assert len(lookup) == len(metadata)
    cohort, cds = [], {}
    for key in ids:
        row = lookup[key]
        assert row['protein_sha256'] == PROTEIN and row['context_sha256'] == CONTEXT and row['longamp'] == 'TRUE'
        raw, code, barcode = row['RNA.sequence'], row['CDS.sequence.longest'], row['Barcode']
        normalized = normalize_rna(raw)
        normalize_rna(code)
        normalize_rna(barcode)
        assert len(raw) == 928 and len(code) == 666 and len(barcode) == 9 and normalized.count(normalize_rna(code)) == 1
        cds[key] = code
        cohort.append({'SequenceID':key, 'split':'development' if key in development else 'holdout',
                       'deposited_full_sha256':text_sha(raw), 'normalized_full_sha256':text_sha(normalized),
                       'CDS_sha256':text_sha(code), 'barcode_sha256':text_sha(barcode),
                       'full_length':len(raw), 'CDS_length':len(code), 'barcode_length':len(barcode)})
    groups = groups_at_distance(development, cds, 10)
    prior_groups = json.loads((BASE/'proposal_review/grouping_metadata.json').read_text())['groups']
    assert groups == prior_groups and len(groups) == 15
    assert sorted(map(len, groups)) == [1]*12+[2,2,6]
    assert RNA.__version__ == '2.7.2' and RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature = 37.0
    assert RNA.params_save(str(HERE/'Turner2004.par'))
    package = Path(RNA.__file__).resolve().parent
    bindings = [Path(RNA.__file__).resolve(), package/'RNA.py']+list(package.glob('_RNA*.pyd'))
    assert len(list(package.glob('_RNA*.pyd'))) == 1
    source_paths = [BASE/name for name in ['permitted_metadata_verified.json','eligibility_results.json',
                    'root_endpoint_proposal.md','report.md','analysis_constructs.tsv','source_manifest.json',
                    'proposal_review/review.md','proposal_review/recomputed_feature_addendum.md',
                    'proposal_review/grouping_metadata.json']]+bindings+[HERE/'Turner2004.par']
    snapshot = {'status':'metadata and model inspection only; no sequence features or cellular targets evaluated',
                'protein_sha256':PROTEIN, 'context_sha256':CONTEXT,
                'all_ids':ids, 'development_ids':development, 'holdout_ids':holdout,
                'development_groups':groups, 'cohort':cohort,
                'metadata_path':str((BASE/'permitted_metadata_verified.json').relative_to(ROOT)),
                'target_table_path':str((BASE/'analysis_constructs.tsv').relative_to(ROOT)),
                'source_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in source_paths},
                'software':{'Python':platform.python_version(),'NumPy':np.__version__,'ViennaRNA':RNA.__version__},
                'vienna':{'version':RNA.__version__,'parameter_set':'Turner2004',
                          'effective_md':snapshot_md(md,RNA), 'pf_float_precision_code':RNA.pf_float_precision(),
                          'precision_interpretation':'0 means double precision, 1 means single precision',
                          'topology':'one linear strand, no cutpoint', 'constraints':'none; no SHAPE or other experimental constraints',
                          'rescaling':'one MFE calculation; one supported exp_params_rescale(mfe) call, sfact=1.07; then one pf call; no manual pf_scale setter, fallback or retry',
                          'probability_tolerance':1e-10,'identity_tolerance':1e-12,
                          'parameter_export_path':str((HERE/'Turner2004.par').relative_to(ROOT)),
                          'unexposed_md_field':'pair_dist is not exposed by this Python binding; retain the pinned package default and do not call comparative-sequence folding',
                          'normalization':'T->U only; keep entire deposited string and experimental barcode; no trimming, cap or poly(A) addition'},
                'inspection_exposure':'Combined source table was hashed as bytes; only its header was inspected. No target or deposited AUP cell was decoded. Metadata grouping and hashes only; no real predictor calculation.'}
    write_new(HERE/'metadata_settings.json',snapshot)
    print(json.dumps({'rows':len(ids),'development':len(development),'holdout':len(holdout),'groups':len(groups),'ViennaRNA':RNA.__version__,'real_folding':False,'outcome_access':False}))


if __name__ == '__main__':
    main()
