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
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3
ITER=("Ala" "Pro" "Pro_Dic")

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

# Loop through Amino Acids
for AA in ${ITER[@]}
do
    # Define input directory
    INPUT=${PROJDIR}/04_rf-norm_${AA}_outs

    # Locate all xml files in the input directory - using find - save to array
    mapfile -t INFILES < <(find "$INPUT" -type f -name '*.xml')

    # Output directory
    FOLD=${PROJDIR}/05_rf-fold_${AA}_outs
    mkdir -p ${FOLD}

    # Loop through each XML file
    for XML in "${INFILES[@]}"
    do
        # Get the directory name from the full path
        dir_path=$(dirname "$XML")
        # Extract the folder name from the directory path
        folder_name=$(basename "$dir_path")
        # Run rf-fold
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
                --output-dir ${FOLD}/${folder_name} \
                ${XML}
    done
done