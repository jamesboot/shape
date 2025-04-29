#!/usr/bin/env bash

#SBATCH --job-name=rf-norm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-6

# Modules
ml purge
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
PAIRS=${PROJDIR}/pairs1.csv

ITER=("Ala" "Pro")

REAGENT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PAIRS} | cut -d ',' -f 1)
TREATED=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PAIRS} | cut -d ',' -f 2)
UNTREATED=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PAIRS} | cut -d ',' -f 3)

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

# Loop through Ala and Pro
for AA in ${ITER[@]}
do
    # Define input directory
    INPUT=${PROJDIR}/03_rf-count_${AA}_outs
    echo "--------------------------------"
    echo "Processing ${AA} samples"
    echo "Input directory: ${INPUT}"
    echo "Reagent: ${REAGENT}"
    echo "Treated: ${TREATED}"
    echo "Untreated: ${UNTREATED}"
    echo "Full path to treated: ${INPUT}/${TREATED}"
    echo "Full path to untreated: ${INPUT}/${UNTREATED}"
    echo "--------------------------------"
    
    # Make output directories
    OUTPUT=${PROJDIR}/04_rf-norm_${AA}_outs
    mkdir -p ${OUTPUT}

    # Check if treated and untreated files exist
    echo "--------------------------------"
    if [[ ! -f "${INPUT}/${TREATED}" ]]; then
        echo "Error: Treated file ${INPUT}/${TREATED} does not exist."
        exit 1
    fi
    if [[ ! -f "${INPUT}/${UNTREATED}" ]]; then
        echo "Error: Untreated file ${INPUT}/${UNTREATED} does not exist."
        exit 1
    fi
    echo "Treated and untreated files exist."
    echo "--------------------------------"

    # If REAGENT is DMS
    if [[ ${REAGENT} == *"DMS"* ]]; then
        # Run rf-norm for DMS
        singularity exec \
            -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
            rf-norm \
                --processors ${SLURM_CPUS_PER_TASK} \
                --overwrite \
                --untreated "${INPUT}/${UNTREATED}" \
                --treated "${INPUT}/${TREATED}" \
                --reactive-bases AC \
                --scoring-method 3 \
                --norm-method 1 \
                --norm-independent \
                --ignore-lower-than-untreated \
                --median-coverage 50 \
                --max-untreated-mut 0.05 \
                --max-mutation-rate 0.2 \
                --output-dir "${OUTPUT}/${TREATED}_vs_${UNTREATED}"
    fi
    
    # If REAGENT is 1M7
    if [[ ${REAGENT} == *"1M7"* ]]; then
        # Run rf-norm for 1M7
        singularity exec \
            -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
            rf-norm \
                --processors ${SLURM_CPUS_PER_TASK} \
                --overwrite \
                --untreated "${INPUT}/${UNTREATED}" \
                --treated "${INPUT}/${TREATED}" \
                --reactive-bases N \
                --scoring-method 3 \
                --norm-method 1 \
                --norm-independent \
                --ignore-lower-than-untreated \
                --median-coverage 50 \
                --max-untreated-mut 0.05 \
                --max-mutation-rate 0.2 \
                --output-dir "${OUTPUT}/${TREATED}_vs_${UNTREATED}"
    fi
done