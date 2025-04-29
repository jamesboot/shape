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
INFILES=(
    "/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape/04_rf-norm_Ala_outs/Ala_DMS.rc_vs_Ala_DMSO.rc/tRNA.xml"
    "/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape/04_rf-norm_Ala_outs/Ala_1M7.rc_vs_Ala_DMSO.rc/tRNA.xml"
    "/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape/04_rf-norm_Pro_outs/ProDic_1M7.rc_vs_ProDic_DMSO.rc/tRNA.xml"
    "/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape/04_rf-norm_Pro_outs/Pro_DMS.rc_vs_Pro_DMSO.rc/tRNA.xml"
    "/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape/04_rf-norm_Pro_outs/Pro_1M7.rc_vs_Pro_DMSO.rc/tRNA.xml"
)
FOLD=${PROJDIR}/05_rf-fold_outs

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

for XML in "${INFILES[@]}"; do
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