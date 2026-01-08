#!/bin/bash
#15-12-2025 Lennart Johansson. Script modified with the help of ChatGPT
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
#chr17_12345678_12345690_15.2 CATCATCATCATCATCATCATCAT	CATCAGTCATCATCATCATCATCAT
#chr18_23456789_23456901_15.2 CATCATCATCATCATCATCATCAT  CATCATCATCATCATCAT
#Here in the first sequence an insertion is included and in the second sequence an interruption is included. In the third sequence shortened allele is created

set -euo pipefail
module load PythonPlus/3.10.4-GCCcore-11.3.0-v23.01.1

############################################
# Configuration
############################################
SCRIPTDIR="/PATH/TO/nanopore/STR_detection/synthetic_controls/scripts"
ANALYSISDIR="/PATH/TO/nanopore/STR_detection/synthetic_controls/v3_UMCGnr_161578_v3_validatie"
CONFIG_TSV="$ANALYSISDIR/insert_synthetic_repeats_config.tsv"
STRAGLRTSV="$ANALYSISDIR/sample01_straglr_1.5.5_VIP_catalog_v1.5.0.tsv"
ORIGINAL_FASTQ="/PATH/TO/nanopore/STR_detection/VIP_8.4.5/sample01/output/intermediates/out_seqtk.fastq.gz"
REFERENCE="/apps/data/vip/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CRAM="/PATH/TO/nanopore/v3_UMCGnr_161578_v3_validatie/run01/results/intermediates/sample01.cram"
CRAI="/PATH/TO/nanopore/v3_UMCGnr_161578_v3_validatie/run01/results/intermediates/sample01.cram.crai"
OUTPUTDIR="$ANALYSISDIR/synthetic_annotate"

mkdir -p "$OUTPUTDIR"

############################################
# Step 1: split FASTQ by allele
############################################
echo "Step 1: splitting fastq by Straglr alleles"

python3 "$SCRIPTDIR/split_fastq_by_straglr_allele_v1.1.py" \
  --outdir "$OUTPUTDIR" \
  --fastq "$ORIGINAL_FASTQ" \
  --tsv "$STRAGLRTSV"

############################################
# Step 2: annotate FASTQ reads
############################################
echo "Step 2: get start and end position and sequence in reads bases on straglr coordinates on alleles defined in config file"
for ALLELE_FASTQ in "$OUTPUTDIR"/*.fastq.gz; do
    python3 "$SCRIPTDIR/annotate_fastq_reads_with_str_coords_v1.4.py" \
      --tsv "$STRAGLRTSV" \
      --fastq "$ALLELE_FASTQ" \
      --cram "$CRAM" \
      --crai "$CRAI" \
      --fasta "$REFERENCE" \
      --outdir "$OUTPUTDIR"
done


############################################
# Step 3: replace STR sequences based on annotation and config TSV
############################################
echo "Step 3: replacing STR sequence with config file new sequence in splitted fastq files"
# Path to FASTQ to modify
FASTQ_STR_ALLELE="$ALLELE_FASTQ"
FASTQ_STR_ALLELE_NAME=$(basename "$FASTQ_STR_ALLELE" .fastq.gz)

# TSV from previous annotation step
ANNOTATION_TSV="$OUTPUTDIR/${FASTQ_STR_ALLELE_NAME}_annotation.tsv"

# Read TSV line by line: base_fastq, original_seq, new_seq
while IFS=$'\t' read -r BASE_FASTQ ORIGINAL_SEQ NEW_SEQ; do
    # skip header or empty lines
    [[ "$BASE_FASTQ" == \#* ]] && continue
    [[ -z "$BASE_FASTQ" ]] && continue

    # filenames based on FASTQ name, not sequence
    ALLELE_FASTQ="$OUTPUTDIR/${BASE_FASTQ}"
    ALLELE_OUT="$OUTPUTDIR/${BASE_FASTQ}.synthetic.fastq.gz"
    LOG_TSV="$OUTPUTDIR/${BASE_FASTQ}.synthetic.log.tsv"

    python3 $SCRIPTDIR/replace_str_sequence_in_fastq_v1.2.py \
      --fastq "$ALLELE_FASTQ" \
      --tsv "$ANNOTATION_TSV" \
      --replace-sequence "$NEW_SEQ" \
      --out-fastq "$ALLELE_OUT" \
      --log-tsv "$LOG_TSV"

done < "$CONFIG_TSV"


############################################
# Step 4: merge synthetic reads back into original FASTQ
############################################
echo "Step 4: Replacing original reads with reads having adaptations in synthetic fastq."
echo "This could take some time, please be patient."
while IFS=$'\t' read -r BASE_FASTQ ORIGINAL_SEQ NEW_SEQ; do
    [[ "$BASE_FASTQ" == \#* ]] && continue
    [[ -z "$BASE_FASTQ" ]] && continue

    SYN_FASTQ="$OUTPUTDIR/${BASE_FASTQ}.synthetic.fastq.gz"
    python3 "$SCRIPTDIR/replace_reads_from_synthetic_fastq.py" \
      --original-fastq "$ORIGINAL_FASTQ" \
      --synthetic-fastq "$SYN_FASTQ" \
      --out-fastq "$OUTPUTDIR/$(basename "$ORIGINAL_FASTQ" .fastq.gz).with_synthetic_${BASE_FASTQ}.fastq.gz" \
      --log-tsv "$OUTPUTDIR/$(basename "$ORIGINAL_FASTQ" .fastq.gz).with_synthetic_${BASE_FASTQ}.simple.log.tsv" \
      --log-detailed-tsv "$OUTPUTDIR/$(basename "$ORIGINAL_FASTQ" .fastq.gz).with_synthetic_${BASE_FASTQ}.log.tsv"
done < "$CONFIG_TSV"
