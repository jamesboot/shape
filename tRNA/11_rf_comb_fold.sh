#!/usr/bin/env bash

# Script to run rf-combine and then rf-fold for tRNA data

#SBATCH --job-name=rf-comb-fold
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='00:10:00'
#SBATCH --mem=32G

# Modules
ml purge
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v4
ITER=("Ala" "Pro" "Pro_Dic")

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

# Variables
AMINO_ACIDS=("Ala" "Pro" "Pro_Dic")
REAGENTS=("DMS" "5NIA")

# Loop through each amino acid
for AA in "${AMINO_ACIDS[@]}"
do

    # Log
    echo "--------------------------------"
    echo "Processing ${AA} samples"
    echo "--------------------------------"
    INPUT=${PROJDIR}/04_rf-norm_${AA}_outs

    # Loop through reagents
    for REAGENT in "${REAGENTS[@]}"
    do
        # Log
        echo "--------------------------------"
        echo "Processing ${AA} with reagent ${REAGENT}"
        echo "--------------------------------"
        # Log
        echo "--------------------------------"
        echo "Making output directories"
        OUTPUT1=${PROJDIR}/06_rf-combine_${AA}_${REAGENT}
        mkdir -p ${OUTPUT1}
        echo "--------------------------------"

        # Run rf-combine
        echo "--------------------------------"
        echo "Running rf-combine for ${AA} with reagent ${REAGENT}"
        echo "--------------------------------"
        singularity exec \
            -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
            rf-combine \
                --processors ${SLURM_CPUS_PER_TASK} \
                --overwrite \
                --stdev \
                --output-dir ${OUTPUT1} \
                ${INPUT}/${REAGENT}_vs_DMSO_0/tRNA.xml ${INPUT}/${REAGENT}_vs_DMSO_1/tRNA.xml ${INPUT}/${REAGENT}_vs_DMSO_2/tRNA.xml

        # Run rf-fold
        echo "--------------------------------"
        echo "Running rf-fold for ${AA} with reagent ${REAGENT}"
        echo "--------------------------------"
        FOLD=${PROJDIR}/07_rf-fold_${AA}_${REAGENT}_comb_outs
        mkdir -p ${FOLD}
        singularity exec \
            -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
            rf-fold \
                --processors ${SLURM_CPUS_PER_TASK} \
                --overwrite \
                --folding-method 1 \
                --img \
                --dotplot \
                --shannon-entropy \
                --windowed \
                --no-lonely-pairs \
                --md 1500 \
                --KT \
                --output-dir ${FOLD} \
                ${OUTPUT1}/tRNA.xml

    done
done    