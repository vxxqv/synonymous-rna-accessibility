"""Selective source reader. Importing this module never opens input files."""
import csv,hashlib,json,math,re,tarfile
from pathlib import Path

IDENTITY=('#CHROM','POS','REF','ALT','POS_NT','REF_NT','ALT_NT','LOCATION','ALT_CODON')
TEXT=('FILTER','UP_REF_NT','DOWN_REF_NT','REF_AA','ALT_AA','AA_POS','MATCHES_MUT_SIG')
COUNTS=('DP','CAO')
NUMERIC=('CAF','R2_MINUS_MED_RP','R1_PLUS_MED_BQ','R2_MINUS_MED_BQ','R2_MINUS_MED_NM')
SELECTED=IDENTITY+TEXT+COUNTS+NUMERIC


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def require_debug():
    if not __debug__:raise RuntimeError('Optimized Python is not permitted')


def checked_json(path, expected):
    if sha(path)!=expected:raise ValueError('Hash mismatch: '+str(path))
    return json.loads(Path(path).read_text(encoding='utf8'))


def parse_summary(stream,expected_header,metadata,contract):
    """Only identities are decoded before the ROI filter; other fields are
    individually allowlisted. Unrequested field bytes are never decoded.
    Every relevant background allele is retained, including nonsynonymous ones.
    """
    first=stream.readline().rstrip(b'\r\n')
    header=first.decode('utf8').split('\t')
    if header!=expected_header or len(header)!=len(set(header)):
        raise ValueError('Unexpected/duplicate source columns')
    index={x:i for i,x in enumerate(header)}
    if not set(SELECTED)<=set(index):raise ValueError('Missing requested columns')
    lo,hi=contract['roi_positions'][str(metadata['ROI'])]
    seen=set();rows=[];excluded=[]
    for number,line in enumerate(stream,2):
        fields=line.rstrip(b'\r\n').split(b'\t')
        if len(fields)!=len(header):raise ValueError('Malformed source row '+str(number))
        get=lambda n:fields[index[n]].decode('utf8')
        chrom=get('#CHROM');position=get('POS')
        if not re.fullmatch(r'\d+',position):raise ValueError('Invalid position')
        pos=int(position)
        if chrom!=contract['chrom'] or not lo<=pos<=hi:
            excluded.append({'line':number,'reason':'outside_fixed_reference_ROI'})
            continue
        row={n:get(n) for n in IDENTITY}
        if row['LOCATION']!='CDS' or not metadata.get('simulation',False) and row['ALT_CODON'] in ('TAA','TGA'):
            excluded.append({'line':number,'reason':'source_LOCATION_or_exact_TAA_TGA_rule'})
            continue
        if not re.fullmatch('[ACGT,]+',row['REF']) or not re.fullmatch('[ACGT,]+',row['ALT']):
            raise ValueError('Unsupported raw allele representation')
        if len(row['REF'])!=len(row['ALT']):raise ValueError('Non-substitution source allele')
        if ',' not in row['REF']:
            if contract['reference'][pos-1:pos-1+len(row['REF'])]!=row['REF']:
                raise ValueError('Raw REF mismatch')
        if not row['POS_NT'].isdigit():raise ValueError('Invalid per-base position')
        posnt=int(row['POS_NT'])
        if row['REF_NT'] not in 'ACGT' or len(row['REF_NT'])!=1 or row['ALT_NT'] not in 'ACGT' or len(row['ALT_NT'])!=1:
            raise ValueError('Invalid per-base allele')
        if posnt<1 or posnt>len(contract['reference']) or contract['reference'][posnt-1]!=row['REF_NT']:
            raise ValueError('Per-base reference mismatch')
        key=(pos,row['REF'],row['ALT'],posnt)
        if key in seen:raise ValueError('Duplicate sample/variant/base row')
        seen.add(key)
        row.update({n:get(n) for n in TEXT})
        for n in COUNTS:
            s=get(n)
            if not re.fullmatch(r'\d+',s):raise ValueError('Invalid nonnegative count '+n)
            row[n]=int(s)
        if row['CAO']>row['DP']:raise ValueError('CAO exceeds DP')
        for n in NUMERIC:
            try:v=float(get(n))
            except ValueError as e:raise ValueError('Missing or malformed required number '+n) from e
            if not math.isfinite(v):raise ValueError('Nonfinite required number '+n)
            if n!='CAF' and v<0:raise ValueError('Negative required classifier feature '+n)
            row[n]=v
        if not 0<=row['CAF']<=1:raise ValueError('CAF outside frequency domain')
        row.update({'POS':pos,'POS_NT':posnt,'Sample':metadata['sample'],
                    'ROI':'ROI_'+str(metadata['ROI']),'source_line':number,
                    'record_id':metadata['sample']+':'+str(number)})
        rows.append(row)
    return rows,excluded


def parse_truth(stream,contract,roi):
    header=None;out=[];seen=set();lo,hi=contract['roi_positions'][str(roi)]
    for number,line in enumerate(stream,1):
        if line.startswith(b'##'):continue
        fields=line.rstrip(b'\r\n').split(b'\t')
        if header is None:
            header=[s.decode('utf8') for s in fields]
            if len(set(header))!=len(header) or not {'#CHROM','POS','REF','ALT','INFO'}<=set(header):raise ValueError('Truth header')
            index={x:i for i,x in enumerate(header)};continue
        if len(fields)!=len(header):raise ValueError('Truth row shape')
        get=lambda n:fields[index[n]].decode('utf8')
        chrom,pos=get('#CHROM'),get('POS')
        if not pos.isdigit():raise ValueError('Truth position')
        if chrom!=contract['chrom'] or not lo<=int(pos)<=hi:continue
        ref,alt,info=get('REF'),get('ALT'),get('INFO')
        first=info.split(';')[0]
        if not re.fullmatch(r'CAO=\d+',first):raise ValueError('Source truth first INFO must be CAO')
        key=(chrom,pos,ref,alt)
        if key in seen:raise ValueError('Duplicate truth identity')
        seen.add(key)
        out.append({'#CHROM':chrom,'POS':int(pos),'REF':ref,'ALT':alt,'INFO':info})
    if header is None:raise ValueError('Missing truth header')
    return out


def write_tsv(path,rows,fields):
    with Path(path).open('x',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def extract(archive,contract,out):
    require_debug();out=Path(out)
    if out.exists():raise FileExistsError('Fresh extraction output required')
    if sha(archive)!=contract['archive_sha256']:raise ValueError('Archive mismatch before source decoding')
    out.mkdir(parents=True)
    allrows=[];receipt=[];active_member=None
    try:
        with tarfile.open(archive,'r:gz') as t:
            names=t.getnames()
            if len(names)!=len(set(names)):raise ValueError('Duplicate archive member')
            allowed=set(contract['members'])|{v[k] for v in contract['simulation'].values() for k in ('summary_member','truth_member')}
            if not allowed<=set(names):raise ValueError('Required source member missing')
            for name in sorted(allowed):
                member=t.getmember(name)
                if not member.isfile() or member.size>30_000_000:raise ValueError('Source member type/size')
            for name in sorted(contract['members']):
                active_member=name;meta=contract['members'][name]
                with t.extractfile(name) as f:rows,exclusions=parse_summary(f,meta['header'],meta,contract)
                allrows.extend(rows);receipt.append({'member':name,'kept_background_rows':len(rows),'identity_exclusions':exclusions})
            fields=list(SELECTED)+['Sample','ROI','source_line','record_id']
            for roi in (1,2):write_tsv(out/f'rows_{roi}.tsv',[r for r in allrows if r['ROI']=='ROI_'+str(roi)],fields)
            for roi,cfg in sorted(contract['simulation'].items()):
                meta={'ROI':int(roi),'sample':f'COMT_{108 if roi=="1" else 109}_train','simulation':True}
                active_member=cfg['summary_member']
                with t.extractfile(active_member) as f:rows,exclusions=parse_summary(f,cfg['header'],meta,contract)
                write_tsv(out/f'simulation_{roi}.tsv',rows,fields)
                active_member=cfg['truth_member']
                with t.extractfile(active_member) as f:truth=parse_truth(f,contract,int(roi))
                write_tsv(out/f'truth_{roi}.tsv',truth,['#CHROM','POS','REF','ALT','INFO'])
                receipt.append({'member':cfg['summary_member'],'simulation_rows':len(rows),'truth_rows':len(truth),'identity_exclusions':exclusions})
    except Exception as error:
        failure={'schema':'COMT_selective_extraction_failure_v1','status':'STOP_NO_FALLBACK',
                 'archive_sha256':contract['archive_sha256'],'active_member':active_member,
                 'completed_members':receipt,'completed_background_rows':len(allrows),
                 'error_type':type(error).__name__,'error':str(error),
                 'exposure_limit':'Required numeric fields may already have been decoded in the active member. Partial-member exposure is not asserted to be zero.'}
        (out/'failure.json').write_text(json.dumps(failure,indent=2)+'\n',encoding='utf8')
        raise
    report={'schema':'COMT_selective_extract_v1','archive_sha256':sha(archive),'background_rows':len(allrows),'members':receipt,'input_universe':'all relevant source-allowed background alleles, independent of125target membership','output_hashes':{p.name:sha(p) for p in sorted(out.glob('*.tsv'))}}
    (out/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
