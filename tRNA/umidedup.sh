#!/usr/bin/env bash

#SBATCH --job-name=umitools_dedup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='6:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-27

# This script performs UMI deduplication on BAM files 

# Load modules
ml purge
ml Anaconda3/2020.07

# Directories
# Edit
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3
DESIGN=${PROJDIR}/samplesheet.csv
# Do not edit
THREADS=${SLURM_CPUS_PER_TASK}
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 1)
ALIGNDIR=${PROJDIR}/02_bowtie2_outs/tRNA_Pro_index
UMIDIR=${PROJDIR}/06_umi_dedup_Pro

# Make dir
mkdir -p ${UMIDIR}

# Extract UMI as the first 10bp of R1.  
# Move extracted sequence to headers and hard-clip the sequence.
echo Starting UMI dedup...
source activate umi_tools_1.1.4
umi_tools dedup \
  --stdin=${ALIGNDIR}/${SAMPLE}_sorted.bam \
  --log=${UMIDIR}/${SAMPLE}.umi_log > ${UMIDIR}/${SAMPLE}_dedup.bam
conda deactivate
echo Finished UMI dedup...
echo "----------------------------------------"