This folder contains the files to create a fastq file with synthetic reads containing adapted repeat sequences.
These can be used in case real positive controls are absent.

The script is made with ONT nanopore data in mind
There are some scripts to be run in preparation:
1. prepare combined fastq file with on-target adaptive sampling file 
run_prepare.sh

2. run vip (scripts not included), to have the cram file.

3. run straglr 1.5.5 with a catalogue for genes of interest to replace with custom sequences
run_straglr_1.5.5_SCAXXX.sh
straglr_1.5.5_benign_DAB1_BEAN1.sh 
VIP_catalog_benign_DAB1_BEAN1.bed

4. run custom scripts to adapt fastq file
input: straglr 1.5.5 output
vip_fam0_FBD12345_straglr_1.5.5_VIP_catalog_v1.5.0.tsv

Scripts for the four steps
annotate_fastq_reads_with_str_coords.py
replace_str_sequence_in_fastq.py
replace_reads_from_synthetic_fastq.py
split_fastq_by_straglr_allele_v1.1.py

Main script to run
insert_synthetic_repeats.sh
Input: tsv file with sequences to replace
insert_synthetic_repeats_config.tsv

#The goal of this script is to use a real fastq sample and replace one or more repeats with synthetic repeats. The inputs are a fastq.gz file and the .tsv outputfile of straglr 1.5.5 (perhaps other
#versions work too) and a configuration file listing which repeats should be replaced with what other repeats. Based on these files the region coordinates are used in step 1 to get the reads per gene and allele into separate fastq files, These will be all alleles
#present in the tsv file,
#This means that when an entire catalog is used, many fastq files are produced, so better is to run straglr, using only the genes that should be replaced.
#This main script combines the four steps, which may make it a bit difficult to create the configuration file (see below), because some knowledge is needed that is easily obtained using intermediate output. However, you can know it beforehand.
#Using the straglr .tsv script, for the genes of interest, look op chrom, start, end and allele (first, second, third and fourteenth (penultimate) column), divided by an underscore, this should make the first column in the configuration file.
#Based on these coordinates, look into a reference genome what is the sequence for that location. This should be entered as the original reference sequence in the second column. In the third column the new alternative sequence should be given.
#This can be an insertion somewhere in the repeat, a created interruption, deletion, or prerepeat or additional repeat.
#If the new sequence contains extra bases the quality of these bases are set to ?, which is equivalent to Q30.

#How to create the configuration file
#insert_synthetic_repeats_config.tsv
# base_fastq_name  original_sequence  new_sequence
#chr16_66490400_66490453_30.4  AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA  AAAATAAAATAAAATAAAATAAGGTAAGGTAAGGTAAGGTAAGGTAAGGTAAGGTAAAATAAAATAAAATAAAATAAAATAAAATAAAA
#chr17_12345678_12345690_15.2 CATCATCATCATCATCATCATCAT  CATCAGTCATCATCATCATCATCAT
#chr18_23456789_23456901_15.2 CATCATCATCATCATCATCATCAT  CATCATCATCATCATCAT
#Here in the first sequence an insertion is included and in the second sequence an interruption is included. In the third sequence shortened allele is created
