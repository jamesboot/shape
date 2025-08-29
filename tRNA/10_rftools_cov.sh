#!/usr/bin/env bash

#SBATCH --job-name=rf-rctools-cov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=ncpu
#SBATCH --time='1:00:00'
#SBATCH --mem=8G

# Modules
ml purge
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3
ITER=("Ala" "Pro")

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

# Loop through Amino Acids
for AA in ${ITER[@]}
do
    # Define input directory
    INPUT=${PROJDIR}/03_rf-count_${AA}_outs

    # Locate all .rc files in the input directory - using find - save to array
    mapfile -t INFILES < <(find "$INPUT" -type f -name '*.rc')

    # Loop through each XML file
    for RC in "${INFILES[@]}"
    do
        # Get the full directory path from the full path
        OUTPATH=$(dirname "$RC")
        # Extract the file name from the OUTPATH
        FILENAME=$(basename "$RC")
        # Run rf-rctools view
        singularity exec \
            -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
            rf-rctools view -t \
                ${RC} > ${OUTPATH}/${FILENAME%.rc}_cov.txt
    done
done