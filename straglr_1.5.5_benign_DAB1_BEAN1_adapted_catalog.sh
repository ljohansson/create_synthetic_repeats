#!/bin/bash

#usage: straglr.py [-h] [--min_ins_size MIN_INS_SIZE] [--min_str_len MIN_STR_LEN] [--max_str_len MAX_STR_LEN] [--min_support MIN_SUPPORT] [--trf_args Match Mismatch Delta PM PI Minscore MaxPeriod] [--chroms CHROMS [CHROMS ...]] [--regions REGIONS] [--exclude EXCLUDE]
#                  [--include_alt_chroms] [--max_cov MAX_COV] [--use_unpaired_clips] [--loci LOCI] [--genotype_in_size] [--include_partials] [--use_mean] [--flank_size FLANK_SIZE] [--sample SAMPLE] [--sex SEX] [--reads_fasta READS_FASTA [READS_FASTA ...]]
#                  [--nprocs NPROCS] [--tmpdir TMPDIR] [--debug] [--version] [--min_cluster_size MIN_CLUSTER_SIZE] [--max_num_clusters MAX_NUM_CLUSTERS] [--max_check_size MAX_CHECK_SIZE] [--min_cluster_d MIN_CLUSTER_D] [--max_bad_cluster_size MAX_BAD_CLUSTER_SIZE]
#                  bam genome_fasta out_prefix

CRAMFOLDER=$1
CRAM="$CRAMFOLDER/$2"
WORKDIR=$3
CATALOGUE="/PATH/TO/resources/VIP_catalog_benign_DAB1_BEAN1_AAAAT_repeat.bed"
REFERENCE="/apps/data/vip/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
CRAMNAME=$(basename ${CRAM%.*})

module purge
module load GCCcore/11.3.0
module load OpenSSL/1.1.1q-GCCcore-11.3.0
module load PythonPlus/3.10.4-GCCcore-11.3.0-v23.01.1
module load zlib/1.3.1
#module load BLAST+/2.16.0-gompi-2024a


if [ ! -d "$WORKDIR/tmp" ]; then
    mkdir "$WORKDIR/tmp"
fi

# Add Python paths
export PYTHONPATH=/PATH/TO/straglr/pysam_dir/lib/python3.10/site-packages:\
/PATH/TO/straglr/pybedtools_310/lib/python3.10/site-packages:\
/PATH/TO/straglr/pathos_310/lib/python3.10/site-packages:$PYTHONPATH

# Add BLAST+ binaries
export PATH=/PATH/TO/straglr/blastn/ncbi-blast-2.17.0+/bin:$PATH



# path to TRF binary
TRF_BIN="/apps/software/TRF/4.09.1-GCCcore-13.3.0/bin/trf"

# call repeats using straglr, ensuring TRF is in PATH
(
  export PATH=$(dirname $TRF_BIN):$PATH

  echo "running straglr"
  python3 /PATH/TO/straglr/straglr-1.5.5/straglr.py \
    $CRAM \
    $REFERENCE \
    ${CRAMNAME}_straglr_1.5.5_VIP_catalog_v1.5.0_adapted_catalog \
    --loci $CATALOGUE \
    --tmpdir $WORKDIR/tmp
)
