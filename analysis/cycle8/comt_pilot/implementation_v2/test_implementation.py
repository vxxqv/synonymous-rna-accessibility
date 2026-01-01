"""Invented fixtures only; never accepts a real archive path."""
import csv,hashlib,io,json,math,random,subprocess,sys,tarfile,tempfile,unittest
from pathlib import Path
from fractions import Fraction as F
import selective_reader as sr
import normalize as nm
import prediction as pr
import driver

P=Path(__file__).resolve().parent


def synthetic_row(sample='S',pos=1,ref='A',alt='C',posnt=None,dp=100,cao=10,line=2):
    return {'#CHROM':'toy','POS':pos,'REF':ref,'ALT':alt,'POS_NT':pos if posnt is None else posnt,
      'REF_NT':'A','ALT_NT':'C','LOCATION':'CDS','ALT_CODON':'CCC','FILTER':'PASS','UP_REF_NT':'A','DOWN_REF_NT':'C',
      'REF_AA':'Q','ALT_AA':'P','AA_POS':'1','MATCHES_MUT_SIG':'True','DP':dp,'CAO':cao,'CAF':cao/dp,
      'R2_MINUS_MED_RP':30,'R1_PLUS_MED_BQ':35,'R2_MINUS_MED_BQ':35,'R2_MINUS_MED_NM':1,
      'Sample':sample,'ROI':'ROI_1','source_line':line,'record_id':sample+':'+str(line)}


def small_contract():
    return {'chrom':'toy','reference':'A'*600,'roi_positions':{'1':[1,300],'2':[301,600]},
      'members':{'S':{'sample':'S','ROI':1,'trait':'R','code_bio':1,'geo_bio':1},
                 'NC':{'sample':'NC','ROI':1,'trait':'NC','code_bio':None,'geo_bio':None}},
      'targets':[{'hgvs_nt':'toy_variant','edit_key':'1:A>C','ROI':1}]}


def bytes_table(rows,poison_unused=True):
    header=list(sr.SELECTED)+['UNREQUESTED']
    lines=['\t'.join(header).encode()]
    for r in rows:
        fields=[str(r[k]).encode() for k in sr.SELECTED]+[b'\xffPOISON' if poison_unused else b'unused']
        lines.append(b'\t'.join(fields))
    return header,io.BytesIO(b'\n'.join(lines)+b'\n')


class Checks(unittest.TestCase):
    def test_selective_poison(self):
        row=synthetic_row();header,f=bytes_table([row]);got,exc=sr.parse_summary(f,header,{'ROI':1,'sample':'S'},small_contract())
        self.assertEqual(got[0]['CAO'],10);self.assertEqual(exc,[])

    def test_outside_identity_does_not_decode_fields(self):
        row=synthetic_row(pos=500);header,f=bytes_table([row]);b=f.getvalue().replace(b'35',b'\xffBAD')
        got,exc=sr.parse_summary(io.BytesIO(b),header,{'ROI':1,'sample':'S'},small_contract())
        self.assertEqual(got,[]);self.assertEqual(len(exc),1)

    def test_duplicate_and_invalid_count(self):
        for rows in ([synthetic_row(),synthetic_row()],[synthetic_row(cao=-1)]):
            header,f=bytes_table(rows)
            with self.assertRaises(ValueError):sr.parse_summary(f,header,{'ROI':1,'sample':'S'},small_contract())

    def test_exact_source_stop_scope(self):
        for alt_codon,simulation,expected in [('TAA',False,0),('TGA',False,0),('TAG',False,1),('TAA',True,1)]:
            row=synthetic_row();row['ALT_CODON']=alt_codon;header,f=bytes_table([row])
            got,_=sr.parse_summary(f,header,{'ROI':1,'sample':'S','simulation':simulation},small_contract());self.assertEqual(len(got),expected)

    def test_classifier_features_nonnegative_domain(self):
        for name in ('R2_MINUS_MED_RP','R1_PLUS_MED_BQ','R2_MINUS_MED_BQ','R2_MINUS_MED_NM'):
            for value in (-1,0,1):
                with self.subTest(feature=name,value=value):
                    row=synthetic_row();row[name]=value;header,stream=bytes_table([row])
                    if value<0:
                        with self.assertRaisesRegex(ValueError,'Negative required classifier feature '+name):
                            sr.parse_summary(stream,header,{'ROI':1,'sample':'S'},small_contract())
                    else:
                        rows,_=sr.parse_summary(stream,header,{'ROI':1,'sample':'S'},small_contract())
                        self.assertEqual(rows[0][name],value)

    def test_normalization_background_and_NC(self):
        rows=[synthetic_row(),synthetic_row(ref='A',alt='G',cao=30,line=3),synthetic_row(sample='NC',cao=0)]
        status={r['record_id']:{'VAR_ID':nm.raw_id(r),'frequency_pass':True,'globally_rejected':False} for r in rows}
        result=nm.normalized(rows,status,small_contract());x=result['measurements'][0]
        expected=F(21,121)-F(1,201)
        self.assertEqual(F(x['corrected_q']),expected)
        self.assertEqual(x['log_q'],math.log(float(expected)))
        status['S:3']['globally_rejected']=True
        result2=nm.normalized(rows,status,small_contract())
        self.assertEqual(F(result2['measurements'][0]['corrected_q']),F(21,181)-F(1,201))

    def test_MNP_median_not_sum(self):
        rows=[synthetic_row(ref='AA',alt='CC',dp=100,cao=10,posnt=1),synthetic_row(ref='AA',alt='CC',dp=102,cao=12,posnt=2,line=3),synthetic_row(ref='A',alt='G',dp=101,cao=30,line=4)]
        c=small_contract();c['targets'][0]['edit_key']='1:A>C;2:A>C'
        st={r['record_id']:{'VAR_ID':nm.raw_id(r),'frequency_pass':True,'globally_rejected':False} for r in rows}
        result=nm.normalized(rows,st,c)
        self.assertEqual(F(result['measurements'][0]['q']),F(23,121));self.assertEqual(result['denominator_lineage'][0]['WT'],'60')

    def test_missing_NC_and_nonpositive(self):
        r=synthetic_row();c=small_contract();st={r['record_id']:{'VAR_ID':nm.raw_id(r),'frequency_pass':True,'globally_rejected':False}}
        z=nm.normalized([r],st,c)['measurements'][0]
        self.assertTrue(z['nc_absence_default_used'])
        nc=synthetic_row(sample='NC',cao=20);st[nc['record_id']]={'VAR_ID':nm.raw_id(nc),'frequency_pass':True,'globally_rejected':False}
        z=nm.normalized([r,nc],st,c)['measurements'][0]
        self.assertIsNone(z['log_q']);self.assertEqual(z['status'],'nonpositive_after_source_NC_subtraction')

    def test_conflicting_DP_and_negative_WT(self):
        for second in [synthetic_row(alt='G',dp=102,cao=10,line=3),synthetic_row(alt='G',dp=100,cao=95,line=3)]:
            rows=[synthetic_row(),second];st={r['record_id']:{'VAR_ID':nm.raw_id(r),'frequency_pass':True,'globally_rejected':False} for r in rows}
            with self.assertRaises(ValueError):nm.normalized(rows,st,small_contract())

    def test_preparation_unknown_not_zero(self):
        c={'members':{},'targets':[{'hgvs_nt':'v','site':'1','ROI':1,'assessment_fold':0,'is_SNV':True}]};measure=[]
        for trait in ('D','R','F','H','M'):
            for b in range(1,5):
                sample=f'{trait}{b}';c['members'][sample]={'sample':sample,'ROI':1,'trait':trait,'code_bio':b,'geo_bio':5-b if trait=='F' else b}
                measure.append({'hgvs_nt':'v','Sample':sample,'log_q':float(b)})
        out=nm.contrast_summaries({'measurements':measure},c)
        self.assertEqual(out['code'][0]['contrasts']['R']['center'],0)
        self.assertEqual(out['code'][0]['contrasts']['R']['lower'],-3)
        measure[0]['log_q']=None
        out=nm.contrast_summaries({'measurements':measure},c)
        self.assertFalse(out['code'][0]['contrasts']['T']['complete']);self.assertTrue(out['code'][0]['contrasts']['L']['complete'])

    def test_prediction_schema_rejects_assessment_target(self):
        k=pr.load_kernel(P.parent/'proxy_models.py')
        with self.assertRaises(ValueError):pr.fit_predict([{'site':'a','R':0,'L':0,'T':0,'ROI':'ROI_1'}],[{'hgvs_nt':'z','site':'z','ROI':'ROI_2','R':0,'L':0,'T':99}],'L',k)

    def test_frozen_CV_exact_load_signal_and_scoring(self):
        k=pr.load_kernel(P.parent/'proxy_models.py');rows=[]
        for i in range(40):
            x=(i%7)-3;z=((i*3)%11)-5;t=2+x+2*z
            cs={n:{'center':float(v),'lower':float(v),'upper':float(v),'complete':True,'resolved_sign':(v>0)-(v<0)} for n,v in [('R',x),('T',t),('L',z),('L_mid',-z)]}
            rows.append({'hgvs_nt':f'toy{i}','site':str(i),'ROI':1 if i%2 else 2,'assessment_fold':i%5,'is_SNV':True,'contrasts':cs})
        data={'code':rows,'geo':rows};pred=pr.cross_validate(data,k);result=pr.evaluate(pred,data,k)
        self.assertEqual(result['cases']['code']['heavy']['scores']['site_MAE']['RNA_load_linear']['float'],0)
        self.assertFalse(result['primary_operational_expansion_pass']) # only40 of125 variants
        altered=json.loads(json.dumps(data));altered['code'][0]['contrasts']['T']['center']=1e6
        pred2=pr.cross_validate(altered,k)
        # All fold0 targets are excluded from every fold0 training argument.
        a={r['hgvs_nt']:r['predictions'] for r in pred['code']['heavy']['predictions'] if r['fold']==0}
        b={r['hgvs_nt']:r['predictions'] for r in pred2['code']['heavy']['predictions'] if r['fold']==0}
        self.assertEqual(a,b)

    def test_optimization_and_missing_authorization(self):
        with tempfile.TemporaryDirectory(dir=P) as td:
            completed=subprocess.run([sys.executable,'-O',str(P/'driver.py'),'--help'],capture_output=True,text=True)
            self.assertNotEqual(completed.returncode,0);self.assertIn('Optimized Python',completed.stderr)
            bad=Path(td)/'auth.json';bad.write_text('{}')
            with self.assertRaises(ValueError):driver.execute('extract',bad,Path(td)/'output')
            self.assertFalse((Path(td)/'output').exists())

    def test_seal_hash_and_overwrite_guards(self):
        with tempfile.TemporaryDirectory(dir=P) as td:
            folder=Path(td);(folder/'predictions.json').write_text('{}');driver.seal(folder,'predict','auth',{})
            driver.stage_manifest(folder/'sealed_manifest.json','predict')
            (folder/'predictions.json').write_text('tampered')
            with self.assertRaises(ValueError):driver.stage_manifest(folder/'sealed_manifest.json','predict')
            with self.assertRaises(FileExistsError):driver.execute('extract',folder/'missing',folder)

    def test_prior_stage_pins_and_unsealed_payload(self):
        with tempfile.TemporaryDirectory(dir=P) as td:
            folder=Path(td);(folder/'predictions.json').write_text('{}')
            pins={name:'pinned_'+name for name in driver.COMMON_PINS}
            driver.seal(folder,'predict','auth',pins)
            driver.stage_manifest(folder/'sealed_manifest.json','predict',pins)
            for name in driver.COMMON_PINS:
                changed={**pins,name:'changed'}
                with self.assertRaises(ValueError):driver.stage_manifest(folder/'sealed_manifest.json','predict',changed)
            bad=json.loads((folder/'sealed_manifest.json').read_text());bad['files']={}
            driver.dump(folder/'sealed_manifest.json',bad)
            with self.assertRaises(ValueError):driver.stage_manifest(folder/'sealed_manifest.json','predict')

    def test_invented_archive_and_failure_exposure_receipt(self):
        with tempfile.TemporaryDirectory(dir=P) as td:
            folder=Path(td);c=small_contract();header,stream=bytes_table([synthetic_row()])
            content=stream.getvalue();c['members']={'S.tsv':{'sample':'S','ROI':1,'header':header}}
            c['simulation']={str(i):{'summary_member':f'sim{i}.tsv','truth_member':f'truth{i}.vcf','header':header} for i in (1,2)}
            contents={'S.tsv':content,'sim1.tsv':content,'sim2.tsv':content,
                      'truth1.vcf':b'#CHROM\tPOS\tREF\tALT\tINFO\ntoy\t1\tA\tC\tCAO=10\n',
                      'truth2.vcf':b'#CHROM\tPOS\tREF\tALT\tINFO\ntoy\t301\tA\tC\tCAO=10\n',
                      'not_allowed.tsv':b'\xffDO_NOT_DECODE'}
            def make_tar(path,parts):
                with tarfile.open(path,'w:gz') as t:
                    for name,data in parts.items():
                        info=tarfile.TarInfo(name);info.size=len(data);t.addfile(info,io.BytesIO(data))
            archive=folder/'invented.tar.gz';make_tar(archive,contents);c['archive_sha256']=sr.sha(archive)
            result=sr.extract(archive,c,folder/'ok');self.assertEqual(result['background_rows'],1)
            self.assertEqual(len(result['members']),3)
            contents['S.tsv']=content.replace(b'\t10\t',b'\tBAD_COUNT\t')
            bad=folder/'bad.tar.gz';make_tar(bad,contents);c['archive_sha256']=sr.sha(bad)
            with self.assertRaises(ValueError):sr.extract(bad,c,folder/'bad_out')
            failure=json.loads((folder/'bad_out/failure.json').read_text())
            self.assertEqual(failure['active_member'],'S.tsv');self.assertEqual(failure['status'],'STOP_NO_FALLBACK')
            self.assertNotIn('receipt.json',{x.name for x in (folder/'bad_out').iterdir()})


if __name__=='__main__':unittest.main()
