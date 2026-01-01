## Step-by-step instructions on processing sequencing reads to raw count

**1. Download sequencing fastq files (BioProject accession no. PRJNA1116243)**
- Example: PTRE-seq library transfected in HeLa cells (biological replicate 1) (accession no. SRR29175799)
   - Contains paired-end fastq files (HELA-1_S23_L001_R1_001.fastq.gz and HELA-1_S23_L001_R2_001.fastq.gz)

**2. Merge paired-end fastq files**
- Command: `bbmerge.sh in1=HELA-1_S23_L001_R2_001.fastq.gz in2=HELA-1_S23_L001_R1_001.fastq.gz out=HELA-1_merged.fastq`
   - Input = HELA-1_S23_L001_R1_001.fastq.gz* and *HELA-1_S23_L001_R2_001.fastq.gz
   - Output = HELA-1_merged.fastq

**3. Find barcodes present in reads and write the results into a demultiplexed text file**
- Command: `python3 demultiplex_barcode.py HELA-1_merged.fastq HELA-1_demux.txt`
   - Input = HELA-1_merged.fastq
   - Output = HELA-1_demux.txt

**4. Split the demultiplexed text file (step 3) into separate report files for each reporter**
- Command: `python3 split_report_by_reporters.py HELA-1_demux.txt`
   - Input = HELA-1_demux.txt & HELA-1_merged.fastq
   - Output = report files for each reporter in `REPORTS/` dir

**5. Split merged fastq files into separate fastq files for each reporter**
- Command: `python3 split_fastq_by_reporter.py HELA-1_demux.txt HELA-1_merged.fastq`
   - Input = HELA-1_demux.txt HELA-1_merged.fastq
   - Output = fastq files for each reporter in `FASTQ_by_reporters/` dir
     
**6. Align separated fastq files (step 5) to corresponding reference sequences**
- Command: `sbatch submit_alignment_bbmap.sh` 
   - Input = separated fastq files in `FASTQ_by_reporters/` dir and reference fasta files in `PTREseq_ref_seq/` dir
   - Output = alignment BAM files for each reporter in `ALIGNMENT/` dir

**7. Categorized reads into full length, spliced, or others based on alignment results**
- Command: `sbatch submit_find_junction_from_alignment.sh`
   - Input = separated report files in `REPORTS/` dir and BAM files in `ALIGNMENT/` dir
   - Output = Updated report files with additional columns indicating the best category for each read

**8. Check for other anomalies (cloning error, mis-priming...)**
- Command: `python3 check_artifacts.py HELA-1`
   - Input = separated report files in `REPORTS/` dir and separated fastq files in `FASTQ_by_reporters/` dir
   - Output = HELA-1_artifacts.txt
     
**9. Concatenate all report text files** 
- Command: `cat REPORTS/HELA-1* > HELA-1_report.txt`

**10. Get raw counts by categories for each reporter**
- Command: `python3 combine_report_and_check_artifacts.py HELA-1`
   - Input = HELA-1_report.txt & HELA-1_artifacts.txt
   - Output = HELA-1_raw_count.txt
- This raw count file can be move to `../ptre_raw_count/` dir. Use `process_raw_count_to_splicing_efficiency.ipynb` notebook to calculate splicing fraction. 


