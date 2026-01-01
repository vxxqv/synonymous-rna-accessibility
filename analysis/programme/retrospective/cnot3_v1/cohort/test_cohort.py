"""Synthetic-only fixtures. Never read assay XLSX or real reference sequences."""
import sys
if not __debug__: raise RuntimeError('Python -O/-OO is refused before input access')
sys.dont_write_bytecode=True
from collections import Counter, defaultdict
from copy import deepcopy
import csv
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('cohort_pipeline_tested',HERE/'pipeline.py')
p=importlib.util.module_from_spec(spec); sys.modules[spec.name]=p; spec.loader.exec_module(p)
c=p.c

def reference(gene,protein='M'):
    aa_codon={aa:min(x for x in c.SENSE if c.CODE[x]==aa) for aa in c.AMINO}
    cds=''.join(aa_codon[aa] for aa in protein)+'TAA'
    cc=Counter(cds[i:i+3] for i in range(0,len(cds)-3,3)); ac=Counter(protein)
    features={**{'codon_'+x:cc[x]/len(protein) for x in c.SENSE},**{'aa_'+a:ac[a]/len(protein) for a in c.AMINO},
              'CDS_length_nt':len(cds),'CDS_GC':sum(x in 'GC' for x in cds)/len(cds),
              'UTR5_length_nt':0,'UTR5_GC':None,'UTR3_length_nt':0,'UTR3_GC':None}
    tid='ENST'+gene[4:]+'.1'
    meta={'stable_gene_id':gene,'transcript_id':tid,'transcript_type':'protein_coding','tags':[],
          'sense_codon_count':len(protein),'CDS_spliced_0based_halfopen':[0,len(cds)]}
    for key,sequence in [('CDS',cds),('transcript',cds),('UTR5',''),('UTR3',''),('protein',protein)]: meta[key+'_sha256']=hashlib.sha256(sequence.encode()).hexdigest()
    return features,meta,cds

def row(gene): return {'GeneID':gene,'HL_WT':2.0,'HL_CNOT3KO':1.0,'LOG2_FC_IP':-3.0,'SE_IP':0.5}

def fixture(n):
    rows=[]; features={}; meta={}; originals=[]
    for i in range(n):
        gene='ENSG'+str(i+1).zfill(11); rows.append(row(gene)); originals.append(gene+'.1')
        features[gene],meta[gene],_=reference(gene)
    return rows,originals,features,meta

def dp(a,b):
    prev=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        now=[i]
        for j,y in enumerate(b,1): now.append(min(prev[j]+1,now[-1]+1,prev[j-1]+(x!=y)))
        prev=now
    return prev[-1]

class Tests(unittest.TestCase):
    def test_strict_original_ids_and_collisions(self):
        ids,report=c.source_membership(['ENSG1.3','ENSG2','ENSG1.3_PAR_Y','ERCC1'])
        self.assertEqual(ids,{'ENSG1','ENSG2'}); self.assertEqual(len(report['unmappable_source_IDs']),2)
        for bad in (['ENSG1.1','ENSG1.2'],['ENSG1','ENSG1'],[123]):
            with self.assertRaises(c.IntegrityError): c.source_membership(bad)

    def test_zero_utrs_negative_changes_retained_no_imputation(self):
        rows,orig,features,meta=fixture(2)
        kept,report=c.eligibility(rows,orig,features,meta)
        self.assertEqual(len(kept),2); self.assertEqual(report['P'],2); self.assertEqual(report['S_over_P'],1)
        self.assertIsNone(features[rows[0]['GeneID']]['UTR5_GC'])
        self.assertFalse(report['UTR_imputation_performed']); self.assertFalse(report['replacement_transcripts_selected'])

    def test_exclusions_no_transcript_replacement_and_sequential_counts(self):
        rows,orig,features,meta=fixture(8); genes=[r['GeneID'] for r in rows]
        orig=orig[1:]; features.pop(genes[1]); meta[genes[2]]['transcript_type']='nonsense_mediated_decay'
        meta[genes[3]]['transcript_type']='protein_coding_LoF'; meta[genes[4]]['tags']=['mRNA_end_NF']
        features[genes[5]]['CDS_GC']=float('nan'); rows[6]['SE_IP']=None; rows[7]['LOG2_FC_IP']=None
        kept,report=c.eligibility(rows,orig,features,meta)
        self.assertEqual(kept,[]); self.assertEqual(len(report['excluded']),8); self.assertEqual(report['sequential_counts'][-1]['remaining'],0)
        self.assertIn('selected_transcript_not_protein_coding',{r['primary_reason'] for r in report['excluded']})
        self.assertEqual(report['S'],7); self.assertEqual(report['A'],0)

    def test_coverage_denominator_and80pct_boundary(self):
        rows,orig,features,meta=fixture(5); features.pop(rows[0]['GeneID'])
        _,a=c.eligibility(rows,orig,features,meta)
        self.assertEqual(a['A_over_S'],.8); self.assertTrue(a['reference_coverage_pass'])
        _,b=c.eligibility(rows,orig[1:],features,meta)
        self.assertEqual(b['S_over_P'],.8); self.assertEqual(b['A_over_S'],1)
        _,empty=c.eligibility([],[],{},{}); self.assertIsNone(empty['A_over_S']); self.assertFalse(empty['reference_coverage_pass'])

    def test_malformed_rows_and_reference_are_integrity_fatal(self):
        rows,orig,f,m=fixture(1)
        for bad in ([rows[0],rows[0]],[{**rows[0],'extra':1}],[{**rows[0],'GeneID':'ENSG1.2'}],[{**rows[0],'HL_CNOT3KO':0}],[{**rows[0],'HL_WT':True}]):
            with self.assertRaises(c.IntegrityError): c.eligibility(bad,orig,f,m)
        broken=deepcopy(f); broken[rows[0]['GeneID']].pop('aa_A')
        with self.assertRaises(c.IntegrityError): c.eligibility(rows,orig,broken,m)

    def test_empty_and_identity_components(self):
        self.assertEqual(c.group_proteins({})['components'],[])
        groups=c.group_proteins({'ENSG3':'M','ENSG1':'M','ENSG2':'A'})
        self.assertEqual(groups['components'],[{'group_ID':'ENSG1','GeneIDs':['ENSG1','ENSG3']},{'group_ID':'ENSG2','GeneIDs':['ENSG2']}])
        self.assertEqual(groups['counts']['exact_pairs_completed'],0)

    def test_transitive_chain_and_boundary_edges(self):
        proteins={'ENSG1':'A'*10,'ENSG2':'C'+'A'*9,'ENSG3':'CC'+'A'*8}
        result=c.group_proteins(proteins)
        self.assertEqual(result['components'],[{'group_ID':'ENSG1','GeneIDs':['ENSG1','ENSG2','ENSG3']}])
        self.assertEqual(len(result['near_protein_edges']),2)
        self.assertEqual(dp(proteins['ENSG1'],proteins['ENSG3']),2)

    def test_multiset_repetition_no_false_negative(self):
        result=c.group_proteins({'ENSG1':'A'*100,'ENSG2':'A'*99})
        self.assertEqual(len(result['components']),1)
        self.assertEqual(result['near_protein_edges'][0]['distance'],1)

    def test_filter_all_pairs_independent_oracle_and_order(self):
        rng=random.Random(7919); sequences=[]
        for _ in range(12):
            base=''.join(rng.choice('ACDE') for _ in range(20)); sequences.append(base)
            sequences.append(('A' if base[0]!='A' else 'C')+base[1:])
            sequences.append(base[:-2])
        proteins={'ENSG'+str(i+1):s for i,s in enumerate(sequences)}
        result=c.group_proteins(proteins)
        reps={}
        for g,s in sorted(proteins.items()): reps.setdefault(s,g)
        expected=set()
        for a,b in itertools.combinations(sorted(reps),2):
            d=dp(a,b)
            if 10*d<=max(len(a),len(b)): expected.add(tuple(sorted([reps[a],reps[b]])))
        actual={tuple(sorted([r['left_representative'],r['right_representative']])) for r in result['near_protein_edges']}
        self.assertEqual(actual,expected)
        reverse=c.group_proteins(dict(reversed(list(proteins.items()))))
        self.assertEqual(result['components'],reverse['components']); self.assertEqual(result['near_protein_edges'],reverse['near_protein_edges'])

    def test_all_resource_caps_fail_closed(self):
        proteins={'ENSG1':'A'*10,'ENSG2':'A'*9+'C'}
        cases=[(c.Limits(0,50000,3_000_000_000,300),'posting_pair_events_preflight'),
               (c.Limits(20_000_000,0,3_000_000_000,300),'candidate_pairs'),
               (c.Limits(20_000_000,50000,99,300),'alignment_cells')]
        for limit,cap in cases:
            with self.assertRaises(c.ResourceLimit) as failure: c.group_proteins(proteins,limits=limit)
            self.assertEqual(failure.exception.receipt['cap'],cap)
            self.assertNotIn('components',failure.exception.receipt); self.assertFalse(failure.exception.receipt['components_complete'])
        exact=c.group_proteins(proteins,limits=c.Limits(20_000_000,1,100,300))
        self.assertEqual(exact['counts']['alignment_cells_completed'],100)

    def test_native_time_overrun_not_accepted_and_checkpoint(self):
        state={'now':0}; progress=[]
        def slow(a,b): state['now']=301; return 1
        with patch.object(c,'exact_distance',slow):
            with self.assertRaises(c.ResourceLimit) as failure:
                c.group_proteins({'ENSG1':'A'*10,'ENSG2':'A'*9+'C'},clock=lambda:state['now'],progress=progress.append)
        self.assertEqual(failure.exception.receipt['cap'],'elapsed_seconds')
        self.assertEqual(failure.exception.receipt['counts']['exact_pairs_completed'],1)
        self.assertTrue(any('pending_alignment_cells' in r['counts'] for r in progress))

    def test_hash_collisions_fail_and_canonical_group_ids(self):
        with self.assertRaises(c.IntegrityError): c.group_proteins({'ENSG1':'M','ENSG2':'A'},digest=lambda s:'same')
        groups=[{'group_ID':'ENSG1','GeneIDs':['ENSG1']},{'group_ID':'ENSG2','GeneIDs':['ENSG2']}]
        with self.assertRaises(c.IntegrityError): c.assign_partitions(groups,digest=lambda s:b'\x01'*32)
        with self.assertRaises(c.IntegrityError): c.assign_partitions([{'group_ID':'ENSG2','GeneIDs':['ENSG1','ENSG2']}])

    def test_full32byte_hash_split_and_heldout_inner_absent(self):
        groups=[{'group_ID':'ENSG'+str(i),'GeneIDs':['ENSG'+str(i)]} for i in range(1,30)]
        records=c.assign_partitions(groups)
        for r in records:
            outer=hashlib.sha256(('cnot3-focus-v1|outer|'+r['group_ID']).encode()).digest()
            self.assertEqual(r['outer'],'holdout' if int.from_bytes(outer,'big')%5==0 else 'training')
            if r['outer']=='holdout': self.assertIsNone(r['inner_fold']); self.assertIsNone(r['inner_SHA256'])
        self.assertEqual(records,c.assign_partitions(list(reversed(groups))))

    def test_split_writer_never_exposes_KO_in_predictors(self):
        rows,orig,features,meta=fixture(2)
        assignment=[{'GeneID':rows[0]['GeneID'],'group_ID':rows[0]['GeneID'],'outer':'training'},
                    {'GeneID':rows[1]['GeneID'],'group_ID':rows[1]['GeneID'],'outer':'holdout'}]
        first=p.split_payloads(rows,features,assignment)
        changed=deepcopy(rows); changed[1]['HL_CNOT3KO']=999
        second=p.split_payloads(changed,features,assignment)
        for key in ('training_predictors','holdout_predictors','training_targets'): self.assertEqual(first[key],second[key])
        self.assertNotEqual(first['holdout_targets'],second['holdout_targets'])
        self.assertNotIn('HL_CNOT3KO',first['holdout_predictors'][0]); self.assertEqual(len(first['holdout_predictors'][0]),92)

    def test_json_collisions_and_authorization_before_row_decode(self):
        with self.assertRaises(c.IntegrityError): p.loads('{"GeneID":"ENSG1","GeneID":"ENSG2"}')
        with tempfile.TemporaryDirectory(dir=HERE) as folder:
            root=Path(folder); authorization=root/'auth.json'
            authorization.write_text('{"real_execution_authorized":false}',encoding='utf-8')
            paths={key:root/'not_opened' for key in ('rows','original_ids','features','metadata','fasta','upstream_gate','protocol')}
            with self.assertRaises(c.IntegrityError): p.execute(paths,authorization,root/'out')
            self.assertFalse((root/'out').exists())

    def test_protein_loader_hashes_duplicates_missing_and_coordinates(self):
        with tempfile.TemporaryDirectory(dir=HERE) as folder:
            root=Path(folder); path=root/'s.fa.gz'; gene='ENSG1'; features,meta,cds=reference(gene,'MA')
            header='>'+meta['transcript_id']+'|'+gene+'.1|x|x|x|x|'+str(len(cds))+'|CDS:1-'+str(len(cds))+'|\n'
            with gzip.open(path,'wt') as stream: stream.write(header+cds+'\n')
            self.assertEqual(p.selected_proteins(path,{gene:meta},[gene],p.sha(path)),{gene:'MA'})
            altered=deepcopy(meta); altered['protein_sha256']='0'*64
            with self.assertRaises(c.IntegrityError): p.selected_proteins(path,{gene:altered},[gene],p.sha(path))
            with self.assertRaises(c.IntegrityError): p.selected_proteins(path,{gene:meta},[gene],'0'*64)
            with gzip.open(path,'wt') as stream: stream.write((header+cds+'\n')*2)
            with self.assertRaises(c.IntegrityError): p.selected_proteins(path,{gene:meta},[gene],p.sha(path))

    def test_5481_gene_end_to_end_synthetic_only(self):
        with tempfile.TemporaryDirectory(dir=HERE) as folder:
            root=Path(folder); paths={key:root/name for key,name in [('rows','r.json'),('original_ids','i.json'),('features','f.csv.gz'),('metadata','m.jsonl.gz'),('fasta','s.fa.gz'),('upstream_gate','g.json'),('protocol','protocol.txt')]}
            rows=[]; originals=[]
            with gzip.open(paths['features'],'wt',newline='') as fs,gzip.open(paths['metadata'],'wt') as ms,gzip.open(paths['fasta'],'wt') as ss:
                writer=csv.DictWriter(fs,fieldnames=['GeneID','transcript_id',*c.FEATURE_NAMES]); writer.writeheader()
                for i,suffix in enumerate(itertools.islice(itertools.product(c.AMINO,repeat=4),5481)):
                    gene='ENSG'+str(i+1).zfill(11); protein='M'+''.join(suffix); values,meta,cds=reference(gene,protein)
                    rows.append(row(gene)); originals.append({'id':gene+'.1'})
                    writer.writerow({'GeneID':gene,'transcript_id':meta['transcript_id'],**values})
                    ms.write(json.dumps(meta)+'\n'); ss.write('>'+meta['transcript_id']+'|'+gene+'.1|x|x|x|x|'+str(len(cds))+'|CDS:1-'+str(len(cds))+'|\n'+cds+'\n')
            paths['rows'].write_text(json.dumps(rows),encoding='utf-8'); paths['original_ids'].write_text(json.dumps({'entries':originals}),encoding='utf-8')
            paths['protocol'].write_text('SYNTHETIC ONLY',encoding='utf-8')
            gate={'schema':'CNOT3_principal_reproduction_gate_v1','status':'passed','cohort_gene_count':5481,'cohort_rows_sha256':p.sha(paths['rows']),
                  'cohort_GeneIDs':sorted(r['GeneID'] for r in rows),'principal_checks':{key:True for key in p.PRINCIPAL_CHECKS},'secondary_warnings':[],
                  'prediction_protocol_sha256':p.sha(paths['protocol'])}
            paths['upstream_gate'].write_text(json.dumps(gate),encoding='utf-8')
            auth={'operation':'CNOT3_cohort_group_partition_v1','real_execution_authorized':True,'reproduction_gate_passed':True,
                  'input_sha256':{key:p.sha(path) for key,path in paths.items()},'code_sha256':{'cohort.py':p.sha(HERE/'cohort.py'),'pipeline.py':p.sha(HERE/'pipeline.py')}}
            authorization=root/'a.json'; authorization.write_text(json.dumps(auth),encoding='utf-8')
            protected=root/'protected_assay.xlsx'; protected.write_bytes(b'POISON: must never be read')
            original_open=Path.open
            def guard(path,*a,**k):
                if str(path).endswith('.xlsx'): raise AssertionError('Protected workbook read')
                return original_open(path,*a,**k)
            with patch.object(p,'REFERENCE_PINS',{key:auth['input_sha256'][key] for key in p.REFERENCE_PINS}),patch.object(Path,'open',guard):
                receipt=p.execute(paths,authorization,root/'o')
            self.assertEqual(receipt['status'],'cohort_group_partition_ready'); self.assertEqual(receipt['final_eligible_genes'],5481)
            guards=p.read_json(root/'o/design_guards.json'); self.assertTrue(all(guards['checks'].values()))
            grouping=p.read_json(root/'o/grouping.json'); self.assertEqual(grouping['counts']['exact_pairs_completed'],0)
            self.assertEqual(len(grouping['components']),5481)
            train=p.read_json(root/'o/model_inputs/training_predictors.json'); heldout=p.read_json(root/'o/model_inputs/holdout_predictors.json')
            self.assertTrue(all('HL_CNOT3KO' not in r and set(r)=={'GeneID','group_ID','HL_WT','LOG2_FC_IP','SE_IP',*c.FEATURE_NAMES} for r in train+heldout))
            self.assertTrue((root/'o/sealed_evaluation/holdout_targets.json').exists())
            self.assertFalse(receipt['training_API_receives_heldout_targets'])

    def test_optimization_flags_fail_before_arguments_or_data(self):
        for script in ('cohort.py','pipeline.py'):
            for flag in ('-O','-OO'):
                result=subprocess.run([sys.executable,flag,'-I','-B',str(HERE/script)],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0); self.assertIn('refused before input access',result.stderr)

if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    receipt={'status':'passed' if result.wasSuccessful() else 'failed','tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
             'scope':'Synthetic fixtures only, including5481 invented genes; no real grouping/partition/assay inputs',
             'code_sha256':{name:p.sha(HERE/name) for name in ('cohort.py','pipeline.py','test_cohort.py')},
             'real_outcome_reads':0,'real_gene_alignments':0,'real_partitions':0}
    with (HERE/'synthetic_receipt.json').open('x',encoding='utf-8') as stream: json.dump(receipt,stream,indent=2); stream.write('\n')
    sys.exit(0 if result.wasSuccessful() else 1)
