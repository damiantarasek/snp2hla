#!/bin/bash

# Basic parameters
max_memory=$1
window_size=$2

# provide a study name
STUDY_NAME=$3

cd /home/mlibydt3/galaxy/tools/

# Set up directory paths - temp output to scratch
TEMP_OUTPUT_DIR=$8
OUTPUT="$TEMP_OUTPUT_DIR"

# Subsetting .... 
DATA=$4
OUT=$4

# Reference panel path.
REF_PANEL=$5

# shortcuts to executables
PLINK=$6
snp2hla=$7

folder_name=`/usr/bin/date +"%Y-%b-%d_%H:%M:%S"`
folder_name=${folder_name}_$3
mkdir -p ${OUTPUT}/${folder_name}

# run the imputation  
$snp2hla ${OUT} ${REF_PANEL} ${OUTPUT}/${STUDY_NAME} $PLINK  $max_memory $window_size

# transpose dosage file
# python myTools/snp2hla/snp2hla_formatter_no_strat_trans.py --dosage ${OUTPUT}/${STUDY_NAME} --out ${OUTPUT}/${STUDY_NAME}.dosage

# transpose dosage file
Rscript $(dirname $0)/transpose_dosage.R ${OUTPUT}/${STUDY_NAME}.dosage ${OUTPUT}/${STUDY_NAME}.fam ${OUTPUT}/${STUDY_NAME}.out.dosage

# copy output to final directory
cp ${OUTPUT}/${STUDY_NAME}* ${OUTPUT}/${folder_name}