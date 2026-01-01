# script to check for experimental and sequencing artifacts

import pandas as pd
import os, sys
from Bio import SeqIO
import edlib

def get_correct_category(category):
	if 'cannot_find_barcode' in category or 'missing_forward_primer' in category or 'RT_error' in category:
		return category
	elif 'full_length' in category:
		return 'full_length'
	elif 'AG_junction_low_mismatch' in category:
		return 'spliced'
	elif 'AG_junction_high_mismatch' in category:
		return 'spliced_ambi'
	elif 'ambiguous' in category:
		return 'ambiguous'
	else:
		return 'chimeric'

# sequence if there was miss priming
query_cannot_find_barcode = 'GGTACCAAACCCGCTG'
query_pum_chimera = 'ATGCATCCGACCGTTGATCTTCCGT'
# sequences if spacer was not cloned
without_spacer = 'ATGCATCGATATCACTCGAG'
with_spacer = 'ATGCATCCGACCGTTGATCTTCCGTGTCAGCTCCGACTACTCGAG'
# sequence if RE modules were not cloned
RE_not_cloned = 'GCTAGCAGATATCTGGTACC'

sample = sys.argv[1]
out = open(f'{sample}_artifacts.txt', 'w')
out.write('read_name\tdist_no_RE\tdist_barcode\tdist_pum_chim\tdist_no_spacer\tdist_spacer\n')

n = 0
for fn in os.listdir('REPORTS/'):
	if fn.startswith(sample):
		df = pd.read_csv(f'REPORTS/{fn}', \
				sep = '\t', \
				usecols = [0, 6], names = ['read_name', 'category'])
		df['best_category'] = df['category'].apply(lambda g: get_correct_category(g))
		
		name = fn.replace('_report.txt', '')
		if not os.path.isfile(f'FASTQ_BY_reporters/{name}_demux.fastq'):
			print(f'{name}_demux.fastq error')
			continue			

		reader = SeqIO.parse(f'FASTQ_BY_reporters/{name}_demux.fastq', 'fastq')
		
		for idx, record in enumerate(reader):
			report_entry = df.iloc[idx]
			read_name = report_entry['read_name']
			if RE_not_cloned in record.seq:
				dist_RE_not_cloned = 0
				dist_barcode = -1
				dist_pumchim = -1
				dist_wo_S = -1
				dist_w_S = -1
			else:	
				if report_entry['best_category'] == 'cannot_find_barcode':
					alignment = edlib.align(query_cannot_find_barcode, record.seq, mode = 'HW', task = 'locations')
					dist_RE_not_cloned = -1
					dist_barcode = alignment['editDistance']
					dist_pumchim = -1
					dist_wo_S = -1
					dist_w_S = -1

				elif report_entry['best_category'] == 'ambiguous':
					dist_RE_not_cloned = -1
					dist_barcode = -1
					alignment = edlib.align(query_pum_chimera, record.seq, mode = 'HW', task = 'locations')
					dist_pumchim = alignment['editDistance']
				
					alignment = edlib.align(without_spacer, record.seq, mode = 'HW')
					dist_wo_S = alignment['editDistance']
					alignment = edlib.align(with_spacer, record.seq, mode = 'HW')
					dist_w_S = alignment['editDistance']

				else:
					dist_RE_not_cloned = -1
					dist_barcode = -1
					dist_pumchim = -1
					dist_wo_S = -1
					dist_w_S = -1

			out.write(f'{read_name}\t{dist_RE_not_cloned}\t{dist_barcode}\t{dist_pumchim}\t{dist_wo_S}\t{dist_w_S}\n')

out.close()
