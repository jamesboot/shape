#!/bin/bash

# This script runs the shapeplotr.R script with specified parameters.

# Input parameters
ENVDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v2
BASEDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape_v4
AMINO_ACIDS=("Ala" "Pro" "Pro_Dic")
REAGENTS=("DMS" "5NIA")
META=${BASEDIR}/meta.csv

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

# Loop through each amino acid
for AA in "${AMINO_ACIDS[@]}"
do

    # Loop through each reagent
    for REAGENT in "${REAGENTS[@]}"
    do
        
        echo "Processing ${AA} with reagent ${REAGENT}"

        # Define the directory containing to fold outputs (for -d and -s options)
        FOLDOUTS=${BASEDIR}/07_rf-fold_${AA}_${REAGENT}_comb_outs
        echo "Fold outputs directory: ${FOLDOUTS}"

        # Define the directory containing the wig files (for -x option)
        WIGOUTS=${BASEDIR}/04_rf-norm_${AA}_outs/${REAGENT}_vs_DMSO
        echo "Wig outputs directory: ${WIGOUTS}"

        # Define the directory containing the rc files
        # Need an exception for Pro_Dic
        if [[ ${AA} == "Pro_Dic" ]]; then
            RCDIR=${BASEDIR}/03_rf-count_Pro_outs
            echo "RC files directory: ${RCDIR}"
        else
            RCDIR=${BASEDIR}/03_rf-count_${AA}_outs
            echo "RC files directory: ${RCDIR}"
        fi
        
        # Define the sample groups for qc based on the meta data
        ARRAYNAME1="${AA}_${REAGENT}"
        ARRAYNAME2="${AA}_DMSO"

        # Use indirect expansion to get the variable names
        eval 'TREATED_GROUP=("${'"$ARRAYNAME1"'[@]}")'
        eval 'UNTREATED_GROUP=("${'"$ARRAYNAME2"'[@]}")'

        # Log
        echo "${AA}_${REAGENT} group samples: ${TREATED_GROUP[@]}"
        echo "${AA}_DMSO group samples: ${UNTREATED_GROUP[@]}"

        # Move to the directory with the r env 
        cd ${ENVDIR}

        # Run the shapeplotr script
        echo "Running shapeplotr for ${AA} with reagent ${REAGENT}"
        ./Rscript /nemo/stp/babs/working/bootj/github/RNA_SHAPE/tRNA/shapeplotr.R \
        -d ${FOLDOUTS}/dotplot/tRNA.dp \
        -x ${WIGOUTS}_0/tRNA.wig,${WIGOUTS}_1/tRNA.wig,${WIGOUTS}_2/tRNA.wig \
        -s ${FOLDOUTS}/shannon/tRNA.wig \
        -o ${BASEDIR}/${AA}_${REAGENT}_shapeplotr.pdf \
        -r "0:75" \
        --limit_y 1.5 \
        --qc \
        --rc_files ${RCDIR}/${TREATED_GROUP[0]}_sorted_cov.txt,${RCDIR}/${TREATED_GROUP[1]}_sorted_cov.txt,${RCDIR}/${TREATED_GROUP[2]}_sorted_cov.txt \
        --rc_controls ${RCDIR}/${UNTREATED_GROUP[0]}_sorted_cov.txt,${RCDIR}/${UNTREATED_GROUP[1]}_sorted_cov.txt,${RCDIR}/${UNTREATED_GROUP[2]}_sorted_cov.txt
        echo "Finished processing ${AA} with reagent ${REAGENT}"
    done
done

