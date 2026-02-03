#!/usr/bin/env bash

#SBATCH --job-name=rf-count
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='1:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml SAMtools/1.18-GCC-12.3.0
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v4
OLDPROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v2

INPUT1=${PROJDIR}/06_umi_dedup_Ala
INPUT2=${PROJDIR}/06_umi_dedup_Pro

INDEXDIR=${OLDPROJDIR}/02_bowtie2_index
INDEX1=${INDEXDIR}/tRNA_Ala_index
FASTA1=${OLDPROJDIR}/tRNA_Ala_AGC_2_1.fa
INDEX2=${INDEXDIR}/tRNA_Pro_index
FASTA2=${OLDPROJDIR}/tRNA_Pro_TGG_3_5.fa

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img
SAMTOOLS=/camp/apps/eb/software/SAMtools/1.18-GCC-12.3.0/bin/samtools

# Run rf-count for tRNA_Ala
OUTDIR1=${PROJDIR}/03_rf-count_Ala_outs
singularity exec \
    -B ${PROJDIR}:${PROJDIR} -B ${INDEXDIR}:${INDEXDIR} -B ${SAMTOOLS}:${SAMTOOLS} ${LEONORE} \
        rf-count \
        --overwrite \
        --fasta ${FASTA1} \
        --count-mutations \
        --min-quality 20 \
        --eval-surrounding \
        --max-collapse-distance 6 \
        --right-deletion \
        --mutation-map \
        --output-dir ${OUTDIR1} \
        ${INPUT1}/*_dedup.bam

# Run rf-count for tRNA_Pro
OUTDIR2=${PROJDIR}/03_rf-count_Pro_outs
singularity exec \
    -B ${PROJDIR}:${PROJDIR} -B ${INDEXDIR}:${INDEXDIR} -B ${SAMTOOLS}:${SAMTOOLS} ${LEONORE} \
        rf-count \
        --overwrite \
        --fasta ${FASTA2} \
        --count-mutations \
        --min-quality 20 \
        --eval-surrounding \
        --max-collapse-distance 6 \
        --right-deletion \
        --mutation-map \
        --output-dir ${OUTDIR2} \
        ${INPUT2}/*_dedup.bam