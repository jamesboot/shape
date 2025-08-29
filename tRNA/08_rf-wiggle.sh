#!/usr/bin/env bash

#SBATCH --job-name=rf-wiggle
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v4
ITER=("Ala" "Pro" "Pro_Dic")

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

# Loop through Amino Acids
for AA in ${ITER[@]}
do
    # Define input directory
    INPUT=${PROJDIR}/04_rf-norm_${AA}_outs
    echo "--------------------------------"
    echo "Processing ${AA} samples"
    echo "Input directory: ${INPUT}"

    # Locate all xml files in the input directory - using find - save to array
    mapfile -t XML_FILES < <(find "$INPUT" -type f -name '*.xml')

    # Log
    echo "Found ${#XML_FILES[@]} XML files in ${INPUT}"
    echo "XML files: ${XML_FILES[@]}"

    # Loop through each XML directory
    for FILE in "${XML_FILES[@]}"
    do
        # Log the current directory being processed
        echo "Processing file: $FILE"
        DIR=$(dirname "$FILE")
        cd "$DIR" # Change to the directory of the XML file
        # Run rf-wiggle for each XML file
        singularity exec \
            -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
            rf-wiggle --overwrite "$FILE"
        echo "Finished processing: $FILE"
        cd "$PROJDIR"  # Return to the project directory
    done
    # Log completion of the current amino acid processing
    echo "Finished processing ${AA} samples"
    echo "--------------------------------"
done