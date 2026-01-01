# script to categorize reads based on alignment results
# input = report text files and sample name (e.g. REPORTS/SH-1_ABAB_2_report.txt and SH-1) 
# output = report text files (see README for description of each column)
 
import os, sys
import pandas as pd
import pysam
import re

fn = sys.argv[1]
sample = sys.argv[2]

def read_fasta(reporter):
    seq = ''
    with open(f'PTREseq_ref_seq/{reporter}.fa', 'r') as fn:
        for line in fn:
             if not line.startswith('>'):
                 seq += line.rstrip()

    return seq

def categorize_read(cigar, dist, ref_seq):

    # covert cigar string into a tuple
    cigar_tuple = re.findall(r'(\d+)([A-Z]|={1})', cigar)
    total_insert = 0
    AG_deletion = 0
    pos_in_ref = 0 # initialize current position

    for length, operation in cigar_tuple:
        # if match, move current position
        if operation == '=' or operation == 'M':
            pos_in_ref += int(length)
        # if deletion, check if it is GT-AG junction
        # then move current position forward
        elif operation == 'D':
            if int(length) > 1:
                start = pos_in_ref
                start_seq = ref_seq[start:start+2]
                end = pos_in_ref + int(length)
                end_seq = ref_seq[end-2:end]
                if start_seq == 'GT' and end_seq == 'AG':
                    AG_deletion += int(length)
                    donor = ref_seq[start-6:start]
                    donor_loc = start 
                    acceptor = ref_seq[end:end+6]
                    acceptor_loc = end
            pos_in_ref += int(length)
        # if insertion, does not move current position forward
        elif operation == 'I' and int(length) > 1:
            total_insert += int(length)
        # ignore soft clipping
        elif operation == 'S':
            pass

    if total_insert > 5:
        # this usually indicate barcode and RE modules mismatched 
        # and is considered PCR chimera
        return 'too_many_inserts' 
    elif AG_deletion > 0:
        # correct edit distance if reads with GT-AG juction
        # for example, if the edit distance is 91, but there is a 90 deletion flanked by GT-AG,
        # then the new edit distance is 91 - 90 = 1
        new_dist = dist - AG_deletion
        junction = f'{donor}_{donor_loc}_{acceptor}_{acceptor_loc}'
        if new_dist < 20:
            # spliced reads
            return f'AG_junction_low_mismatch-{junction}'
        else:
            # spliced reads, but ambiguous
            return f'AG_junction_high_mismatch-{junction}'
    else:
        if dist < 20:
            # full length reads
            return 'full_length'
        else:
            # reads with ambiguous deletions
            return 'ambiguous'

def apply_categorize_read(read_name, read_type, reporters, _cigar, _dist, ref_seq):
    if read_type == 'perfect_match' or read_type == 'one_mismatch':
        if _cigar is None or _dist is None:
            category = 'one_more_errors' 
        else:    
            category = categorize_read(_cigar, _dist, ref_seq)

        cigar = _cigar
        dist = _dist
    else:
        cigar = read_type
        dist = read_type
        category = read_type
    return cigar, dist, category

linestowrite = []
with open(f'REPORTS/{fn}', 'r') as report:
	if 'error' in fn:
		# for the error report, add the error to the 4 new columns
		for line in report:
			spl = line.rstrip().split('\t')
			original_line = '\t'.join(spl[:4])
			error = spl[2]
			towrite = f'{original_line.rstrip()}\t{error}\t{error}\t{error}\t{error}\n'
			linestowrite.append(towrite)	
	else:
        # get reporter name
		reporter = fn.split('_')[1] + '_' + fn.split('_')[2]
		# read alignment BAM file
		bam_file = f'ALIGNMENT/{sample}_{reporter}_demux.sam.bam'
		bam = pysam.AlignmentFile(bam_file, 'r')
		bam_dict = {}

        # key = read ID; value = cigar string and tag
		for query in bam.fetch(until_eof = True):
			name = query.query_name.split(' ')[0] 
			cigar = query.cigarstring
			try:
				tag = query.get_tag('NM')
			except:
				tag = None
			bam_dict[name] = (cigar, tag)	
			

		# read reference sequence
		ref_seq = read_fasta(reporter)    

        # read each line of the report and add 4 new columns: cigar string, edit distance, category, and reporter
		for line in report:	
			spl = line.rstrip().split('\t')
			original_line = '\t'.join(spl[:4])
			_cigar, _dist = bam_dict[spl[0]]
			if _cigar is None or _dist is None:
				cigars, distances, category = ['error', 'error', 'error']
			else:	
				cigars, distances, category = apply_categorize_read(spl[0], spl[2], spl[3], _cigar, _dist, ref_seq)
			
			towrite = f'{original_line.rstrip()}\t{cigars}\t{distances}\t{category}\t{reporter}\n'	
			linestowrite.append(towrite)	

with open(f'REPORTS/{fn}', 'w') as out:
	for line in linestowrite:
		out.write(line)

