#!/usr/bin/env bash

#SBATCH --job-name=rf-count
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml SAMtools/1.18-GCC-12.3.0
ml Bowtie2/2.5.1-GCC-12.3.0
ml RNAFramework/2.9.0-foss-2021b

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
FASTA=${PROJDIR}/T7_tRNA_Ala_AGC_2_1.fasta
OUTDIR1=${PROJDIR}/03_rf-map_outs
OUTDIR2=${PROJDIR}/03_rf-count_outs
INPUTDIR=${PROJDIR}/01_preprocess_reads_outs/05_adjusted_header
INDEX=${PROJDIR}/02_bowtie2_index/tRNA_index

# Run rf-map
# cq5 = 5'-end quality trimming must be avoided when analyzing data from RT-stop-based methods?
rf-map \
    -b2 \
    -cqo \
    -cq5 20 \
    -bs \
    -bl 15 \
    -bN 1 \
    -bD 20 \
    -bR 3 \
    -bdp 100 \
    -bma 2 \
    -bmp 6,2 \
    -bdg 5,1 \
    -bfg 5,1 \
    -mp "--maxins 200" \
    -bi ${INDEX} \
    -o ${OUTDIR1} \
    -ow \
    ${INPUTDIR}/*.fq.gz

# Run rf-count
rf-count \
    --fast \
    -o ${OUTDIR2} \
    -ow \
    -f ${FASTA} \
    -m \
    -na \
    -md 10 \
    ${OUTDIR1}/*.bam



