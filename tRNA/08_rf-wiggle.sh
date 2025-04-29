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
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
ITER=("Ala" "Pro")

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

# Loop through Ala and Pro
for AA in ${ITER[@]}
do
    # Define input directory
    INPUT=${PROJDIR}/04_rf-norm_${AA}_outs
    echo "--------------------------------"
    echo "Processing ${AA} samples"
    echo "Input directory: ${INPUT}"
    echo "--------------------------------"

    # Create .wig files
    find "${INPUT}" -type f -name "*.xml" | while read -r XML_FILE; do
    WIG_FILE="${XML_FILE%.xml}.wig"
    if [[ -f "$WIG_FILE" ]]; then
        echo "Skipping wiggle for $XML_FILE — .wig exists."
        continue
    fi

    echo "Generating wiggle: $XML_FILE"
    XML_DIR=$(dirname "$XML_FILE")
    cd "$XML_DIR" || exit

    singularity exec \
        -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-wiggle --overwrite "$XML_FILE"

    cd "$PROJECT" || exit
    done
done