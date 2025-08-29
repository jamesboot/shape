#!/usr/bin/env bash

#SBATCH --job-name=rf-merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml SAMtools/1.18-GCC-12.3.0
ml Singularity/3.6.4

# Input parameters
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v3
INPUT1=${PROJDIR}/03_rf-count_Pro_outs

# Static directories
LEONORE=/camp/lab/bauerd/home/shared/singularity/amchakrabarti-leonore-0.1.0.img

G1=(${INPUT1}/SAN6476A99_sorted.rc ${INPUT1}/SAN6476A100_sorted.rc ${INPUT1}/SAN6476A101_sorted.rc)
G1NAME=Ala_DMS
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G1[@]} -o ${INPUT1}/${G1NAME}.rc

G2=(${INPUT1}/SAN6476A102_sorted.rc  ${INPUT1}/SAN6476A103_sorted.rc ${INPUT1}/SAN6476A104_sorted.rc)
G2NAME=Pro_DMS
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G2[@]} -o ${INPUT1}/${G2NAME}.rc

G3=(${INPUT1}/SAN6476A105_sorted.rc ${INPUT1}/SAN6476A106_sorted.rc ${INPUT1}/SAN6476A107_sorted.rc)
G3NAME=ProDic_DMS
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G3[@]} -o ${INPUT1}/${G3NAME}.rc

G4=(${INPUT1}/SAN6476A108_sorted.rc ${INPUT1}/SAN6476A109_sorted.rc ${INPUT1}/SAN6476A110_sorted.rc)
G4NAME=Ala_1M7
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G4[@]} -o ${INPUT1}/${G4NAME}.rc

G5=(${INPUT1}/SAN6476A111_sorted.rc ${INPUT1}/SAN6476A112_sorted.rc ${INPUT1}/SAN6476A113_sorted.rc)
G5NAME=Pro_1M7
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G5[@]} -o ${INPUT1}/${G5NAME}.rc

G6=(${INPUT1}/SAN6476A114_sorted.rc ${INPUT1}/SAN6476A115_sorted.rc ${INPUT1}/SAN6476A116_sorted.rc)
G6NAME=ProDic_1M7
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G6[@]} -o ${INPUT1}/${G6NAME}.rc

G7=(${INPUT1}/SAN6476A117_sorted.rc ${INPUT1}/SAN6476A118_sorted.rc ${INPUT1}/SAN6476A119_sorted.rc)
G7NAME=Ala_DMSO
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G7[@]} -o ${INPUT1}/${G7NAME}.rc

G8=(${INPUT1}/SAN6476A120_sorted.rc ${INPUT1}/SAN6476A121_sorted.rc ${INPUT1}/SAN6476A122_sorted.rc)
G8NAME=Pro_DMSO
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G8[@]} -o ${INPUT1}/${G8NAME}.rc

G9=(${INPUT1}/SAN6476A123_sorted.rc ${INPUT1}/SAN6476A124_sorted.rc ${INPUT1}/SAN6476A125_sorted.rc)
G9NAME=ProDic_DMSO
singularity exec \
    -B ${PROJDIR}:${PROJDIR} ${LEONORE} \
        rf-rctools merge ${G9[@]} -o ${INPUT1}/${G9NAME}.rc