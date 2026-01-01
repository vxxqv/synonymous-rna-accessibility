"""Verify the published symbolic design, never read source phenotype data."""
from pathlib import Path
from itertools import product
from math import prod
from Bio.Seq import Seq
import json,hashlib

if not __debug__:raise RuntimeError('Optimized Python prohibited')
p=Path('work/research_cycle8/neighborhood_gate')
dest=p/'structure_and_provenance.json'
if dest.exists():raise FileExistsError(dest)
patterns=['TCN','AAR','GGN','GAR','GAR','CTN','TTY','ACN']
iupac={'N':'ACGT','R':'AG','Y':'CT'}
domains=[[''.join(w) for w in product(*(iupac.get(c,c) for c in codon))] for codon in patterns]
assert [len(x) for x in domains]==[4,2,4,2,2,4,2,4]
translations={str(Seq(''.join(codons)).translate(table=1)) for codons in product(*domains)}
assert translations=={'SKGEELFT'}
assert all(sum(a!=b for a,b in zip(x,y))==1 for d in domains for x in d for y in d if x!=y)
n=prod(map(len,domains));degree=sum(len(d)-1 for d in domains)
assert (n,degree,n*degree//2)==(4096,16,32768)
manifest={'scope':'Symbolic intended design only; no measured genotype/phenotype graph certified',
          'schmitz':{'published_variable_codons_after_ATG':patterns,'codon_domains':domains,
                      'eight_residue_translation':sorted(translations),'verified_variable_sequences':n,
                      'graph':'K4^4 x K2^4','SNV_degree':degree,'intended_edges':n*degree//2,
                      'measured_complete_neighborhoods':None,'deposited_variant_values_read':False,
                      'complete_construct_sequence_authenticated':False,'analysis_code_version_authenticated':False,
                      'primary_share':'https://figshare.com/s/51c4820b2bee85c94007',
                      'previous_share_receipt':'work/joint_function_cycle5/splicing_sources/retrieval_initial.json',
                      'previous_access_status':'202 empty; not retried in this task'},
          'tp53':{'source_DOI':'10.1158/1541-7786.MCR-17-0245','raw_accession':'PRJNA384242',
                  'data_file_map':{'Supplementary_File_1':'NIHMS888465-supplement-7.gz',
                                   'Supplementary_File_2_notebook':'NIHMS888465-supplement-8.txt',
                                   'Supplementary_Table_1_oligonucleotides':'NIHMS888465-supplement-5.xlsx'},
                  'attachment_access':'403; preserved and not retried','processed_schema_authenticated':False,
                  'multi_mutant_retention_authenticated':False,'measured_complete_neighborhoods':None,
                  'analysis_notebook_version_authenticated':False,'source_phenotype_values_read':False,
                  'full_native_protein_assay':False,'reporter':'TP53 residues97-292 fused to EGFP with exon5-7 processing context'},
          'history':'Fixed seven-source F0 audit, two sources owned here; not exhaustive prevalence evidence',
          'exposure':'Primary Methods and already published procedure/coverage/filter/caption summaries; no source outcome cells',
          'acquisition_failures':['PMC browser challenges, no interaction','TP53 EuropePMC fullTextXML500',
                                  'TP53 publisher403','Indexed author PDF404','Three explicit PMC media403 responses'],
          'no_new_outcome_analysis':True,'historical_git_reauthenticated':False,
          'files':[]}
for f in sorted(p.iterdir()):
    if f.is_file() and f.name!=dest.name:
        manifest['files'].append({'path':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),
                                   'local_source_cache_not_automatic_archive_whitelist':f.suffix in ['.xml','.docx','.png','.xlsx','.gz','.txt']})
assert '\u2014' not in (p/'report.md').read_text(encoding='utf-8')
dest.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'report_sha256':hashlib.sha256((p/'report.md').read_bytes()).hexdigest(),
                  'manifest_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
                  'status':'F0_TWO_SOURCE_GATE_COMPLETE_NO_OUTCOMES'}))
