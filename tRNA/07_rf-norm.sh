#!/usr/bin/env bash

#SBATCH --job-name=rf-norm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3
META=${PROJDIR}/meta.csv

ITER=("Ala" "Pro")

# Define sample groups
# DMS samples
mapfile -t Ala_DMS < <(grep -w 'Ala' $META | grep -w 'DMS' | cut -d ',' -f 1)
mapfile -t Pro_DMS < <(grep -w 'Pro' $META | grep -w 'DMS' | cut -d ',' -f 1)
mapfile -t Pro_Dic_DMS < <(grep -w 'Pro_Dic' $META | grep -w 'DMS' | cut -d ',' -f 1)
# 5NIA samples
mapfile -t Ala_5NIA < <(grep -w 'Ala' $META | grep -w '5NIA' | cut -d ',' -f 1)
mapfile -t Pro_5NIA < <(grep -w 'Pro' $META | grep -w '5NIA' | cut -d ',' -f 1)
mapfile -t Pro_Dic_5NIA < <(grep -w 'Pro_Dic' $META | grep -w '5NIA' | cut -d ',' -f 1)
# DMSO samples
mapfile -t Ala_DMSO < <(grep -w 'Ala' $META | grep -w 'DMSO' | cut -d ',' -f 1)
mapfile -t Pro_DMSO < <(grep -w 'Pro' $META | grep -w 'DMSO' | cut -d ',' -f 1)
mapfile -t Pro_Dic_DMSO < <(grep -w 'Pro_Dic' $META | grep -w 'DMSO' | cut -d ',' -f 1)

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
    # Need an exception for Pro_Dic
    if [[ ${AA} == "Pro_Dic" ]]; then
        # Define input directory
        INPUT=${PROJDIR}/03_rf-count_Pro_outs
        # Make output directories
        OUTPUT=${PROJDIR}/04_rf-norm_${AA}_outs
        mkdir -p ${OUTPUT}
    else
        # Define input directory
        INPUT=${PROJDIR}/03_rf-count_${AA}_outs
        # Make output directories
        OUTPUT=${PROJDIR}/04_rf-norm_${AA}_outs
        mkdir -p ${OUTPUT}
    fi

    # Loop through reagents
    for REAGENT in "${REAGENTS[@]}"
    do

        # Log
        echo "--------------------------------"
        echo "Processing ${AA} with reagent ${REAGENT}"
        # Define treated and untreated samples based on amino acid and reagent
        ARRAYNAME1="${AA}_${REAGENT}"
        ARRAYNAME2="${AA}_DMSO"
        # Use indirect expansion to get the variable names
        eval 'TREATED_GROUP=("${'"$ARRAYNAME1"'[@]}")'
        eval 'UNTREATED_GROUP=("${'"$ARRAYNAME2"'[@]}")'

        # Log
        echo "Treated group samples: ${TREATED_GROUP[@]}"
        echo "Untreated group samples: ${UNTREATED_GROUP[@]}"
        echo "--------------------------------"

        # Get the length of the arrays
        LENGTH1=${#TREATED_GROUP[@]}
        LENGTH2=${#UNTREATED_GROUP[@]}
        # Check if both arrays have the same length
        if [[ ${LENGTH1} -ne ${LENGTH2} ]]; then
            echo "Error: Treated and untreated groups have different lengths."
            exit 1
        fi

        # Loop through numbers 1 to length of array
        for (( NUM=0; NUM<LENGTH1; NUM++ ))
        do

            # Log
            echo "--------------------------------"
            TREATED_SAMPLE=${TREATED_GROUP[$NUM]}
            UNTREATED_SAMPLE=${UNTREATED_GROUP[$NUM]}
            echo "Replicate number: ${NUM}"
            echo "Processing sample: ${TREATED_SAMPLE} vs ${UNTREATED_SAMPLE}"
            echo "Input directory: ${INPUT}"
            echo "Full path to treated: ${INPUT}/${TREATED_SAMPLE}_L008_sorted.rc"
            echo "Full path to untreated: ${INPUT}/${UNTREATED_SAMPLE}_L008_sorted.rc"
            echo "--------------------------------"

            # Check if treated and untreated files exist
            echo "--------------------------------"
            if [[ ! -f "${INPUT}/${TREATED_SAMPLE}_L008_sorted.rc" ]]; then
                echo "Error: Treated file ${INPUT}/${TREATED_SAMPLE}_L008_sorted.rc does not exist."
                exit 1
            fi
            if [[ ! -f "${INPUT}/${UNTREATED_SAMPLE}_L008_sorted.rc" ]]; then
                echo "Error: Untreated file ${INPUT}/${UNTREATED_SAMPLE}_L008_sorted.rc does not exist."
                exit 1
            fi
            echo "Treated and untreated files exist."
            echo "--------------------------------"

            # If REAGENT is DMS
            if [[ ${REAGENT} == *"DMS"* ]]; then

                echo "--------------------------------"
                echo "Running rf-norm for DMS"
                # Run rf-norm for DMS
                singularity exec \
                    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
                    rf-norm \
                        --processors ${SLURM_CPUS_PER_TASK} \
                        --overwrite \
                        --untreated "${INPUT}/${UNTREATED_SAMPLE}_L008_sorted.rc" \
                        --treated "${INPUT}/${TREATED_SAMPLE}_L008_sorted.rc" \
                        --reactive-bases AC \
                        --scoring-method 3 \
                        --norm-method 1 \
                        --norm-independent \
                        --ignore-lower-than-untreated \
                        --median-coverage 50 \
                        --max-untreated-mut 0.05 \
                        --max-mutation-rate 0.2 \
                        --output-dir "${OUTPUT}/${REAGENT}_vs_DMSO_${NUM}"
                echo "rf-norm for DMS completed."
                echo "--------------------------------"

            fi
    
            # If REAGENT is 5NIA
            if [[ ${REAGENT} == *"5NIA"* ]]; then

                echo "--------------------------------"
                echo "Running rf-norm for 5NIA"
                # Run rf-norm for 5NIA
                singularity exec \
                    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
                    rf-norm \
                        --processors ${SLURM_CPUS_PER_TASK} \
                        --overwrite \
                        --untreated "${INPUT}/${UNTREATED_SAMPLE}_L008_sorted.rc" \
                        --treated "${INPUT}/${TREATED_SAMPLE}_L008_sorted.rc" \
                        --reactive-bases N \
                        --scoring-method 3 \
                        --norm-method 1 \
                        --norm-independent \
                        --ignore-lower-than-untreated \
                        --median-coverage 50 \
                        --max-untreated-mut 0.05 \
                        --max-mutation-rate 0.2 \
                        --output-dir "${OUTPUT}/${REAGENT}_vs_DMSO_${NUM}"
                echo "rf-norm for 5NIA completed."
                echo "--------------------------------"

            fi

        done # End of NUM loop

    done # End of REAGENT loop

done # End of AA loop