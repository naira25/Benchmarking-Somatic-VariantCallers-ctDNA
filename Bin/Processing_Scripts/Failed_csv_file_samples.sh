#!/bin/bash

# Create a Sample Sheet with fastq pairs for each sample in csv format
# Comma separated values (csv) file with the minimal configuration for the samples sheet as; Patient, Sample, Lane, fastq_1, fastq_2
# The characteristics are; at least 3 columns and header line

samples=("B" "D" "E")

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA/"
OUTPUT_DIR="$BASE_DIR/Metadata"

CSV_FILE="$OUTPUT_DIR/samplesheet.csv"

echo "patient,status,sample,lane,fastq_1,fastq_2" > $CSV_FILE

for sample in ${samples[@]}; do
    echo "Processing Sample ${sample}"

    SAMPLE_DIR=$BASE_DIR/Test_Samples/Sample_${sample}/fastq_data/fastq
    
    for file in $SAMPLE_DIR/*_1.fastq.gz; do
        filename=$(basename "$file")

        Patient="Sample_${sample}"
        
        if [[ ${sample} == *"B"* ]]; then
            Status="0"
        else
            Status="1"
        fi

        Sample_ID=${filename%_1.fastq.gz} 
        Lane="1"
        fastq_1=$(realpath "$file")
        fastq_2=$(realpath "${file/_1.fastq.gz/_2.fastq.gz}")

    echo "$Patient,$Status,$Sample_ID,$Lane,$fastq_1,$fastq_2" >> $CSV_FILE
    done
done
