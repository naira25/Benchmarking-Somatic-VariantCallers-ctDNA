#!/bin/bash

samples=("B" "D" "E") # Define an array of sample names

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
INPUT_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA/Samples" 
OUTPUT_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA/Samples"
NEXTFLOW_BIN="$BASE_DIR/nextflow"

# Loop through each sample in the array
for sample in "${samples[@]}"; do
    echo "Processing Sample ${sample}"
    
    SAMPLE_OUT="$OUTPUT_DIR/Sample_${sample}"
    mkdir -p "$SAMPLE_OUT"
    
    # Definim on es guardarà el fitxer d'IDs (dins de la carpeta de la mostra)
    INPUT_TXT="$SAMPLE_OUT/SRR_Acc_${sample}.txt"

    {
        echo "ids"
        awk -F',' '/_ROC2_/ && /(ST10|ST20|ST21)/ && /25ng/ {gsub(/"/, "", $1); print $1}' "${INPUT_DIR}/SRR_Acc_${sample}.csv"
    } > "$INPUT_TXT"

    # Run the Nextflow pipeline for each sample using the generated text file as input
    $NEXTFLOW_BIN run nf-core/fetchngs \
        --input "$INPUT_TXT" \
        --download_method aspera \
        --outdir "$SAMPLE_OUT/fastq_data" \
        -profile docker \
        -work-dir "$SAMPLE_OUT/work" \
        -resume
done
