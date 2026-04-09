#!/bin/bash

subsamples=("1" "2" "3" "4" "5" "6")

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA/"
OUTPUT_DIR="$BASE_DIR/Metadata/SampleSheets_Subsamples"

mkdir -p "$OUTPUT_DIR"

for sub in ${subsamples[@]}; do
    echo "Processing Subsample ${sub}"

    # CORRECCIÓ: Fem servir ${sub} per tenir un fitxer per cada carpeta
    CSV_FILE="$OUTPUT_DIR/samplesheet_${sub}.csv"
    
    echo "patient,status,sample,lane,fastq_1,fastq_2" > $CSV_FILE

    SAMPLE_DIR="$BASE_DIR/Subsamples/Subsample_${sub}"
    
    for file in $SAMPLE_DIR/*_1.fastq.gz; do

        filename=$(basename "$file")

        Patient="Patient_${sub}"
        
        if [ -f "$BASE_DIR/Data/Samples/Sample_B/fastq_data/fastq/$filename" ]; then
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