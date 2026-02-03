#!/usr/bin/env bash

# A bash script to submit via sbatch to a slurm scheduler that merges R1 and R2 reads of the same sample from different lanes that exist in the same directory.

#SBATCH --job-name=merge_lanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Directories
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3
FASTQDIR=/nemo/stp/sequencing/inputs/instruments/fastq/20250722_LH00442_0158_A232WTHLT3/fastq/PM23265
OUTDIR=${PROJDIR}/merged_reads

# Create output directory if it does not exist
mkdir -p ${OUTDIR}

# Merge R1 and R2 reads from different lanes for each sample
for SAMPLE in $(ls ${FASTQDIR} | grep -oP '^[^_]+' | sort -u); do
    # Find all R1 and R2 files for the sample
    R1_FILES=$(ls ${FASTQDIR}/${SAMPLE}*_R1*.fastq.gz)
    R2_FILES=$(ls ${FASTQDIR}/${SAMPLE}*_R2*.fastq.gz)
    # Merge R1 files
    if [ -n "${R1_FILES}" ]; then
        cat ${R1_FILES} > ${OUTDIR}/${SAMPLE}_R1_001.fastq.gz
        echo "Merged R1 reads for sample ${SAMPLE} into ${OUTDIR}/${SAMPLE}_R1_001.fastq.gz"
    else
        echo "No R1 files found for sample ${SAMPLE}"
    fi
    # Merge R2 files
    if [ -n "${R2_FILES}" ]; then
        cat ${R2_FILES} > ${OUTDIR}/${SAMPLE}_R2_001.fastq.gz
        echo "Merged R2 reads for sample ${SAMPLE} into ${OUTDIR}/${SAMPLE}_R2_001.fastq.gz"
    else
        echo "No R2 files found for sample ${SAMPLE}"
    fi
done
