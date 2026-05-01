#!/bin/bash

# ==================================================================================================================================
# FASTQ DATA RETRIEVAL
# ==================================================================================================================================
# This script retrieves raw sequencing data (FASTQ files) for Samples B, D, and E
# with the nf-core/fetchngs Nextflow pipeline using the following steps:
# 1. Preparing the input metadata with SRA accession IDs.
# 2. Downloading raw reads from public repositories based on the provided accessions
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA/"
INPUT_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA/samples"
NEXTFLOW_BIN="$BASE_DIR/nextflow"

samples=("B" "D" "E") # Define an array with the samples to process

# Iterate though each sample to fecth FASTQ files
for sample in "${samples[@]}"; do

    echo "Processing Sample ${sample}"
    
    # Create the sample-specific output directories for the downloaded FASTQ files
    SAMPLE_DIR="$INPUT_DIR/Sample_${sample}"
    mkdir -p "$SAMPLE_DIR"
    
    # 1. Convert the .txt file with SRA accession IDs to a .csv
    INPUT_TXT="$INPUT_DIR/SRR_Acc_${sample}.txt" # Define the path to the input .txt file 
    INPUT_CSV="$INPUT_DIR/SRR_Acc_${sample}.csv" # Define the path to the input .csv file
    cat "$INPUT_TXT" > "$INPUT_CSV"

    # 2. Run the nf-core/fetchngs Nextflow pipeline with the following options:
    # --input: Provide the input .csv file with SRA accession IDs
    # --outdit: Define the output directory for the downloaded FASTQ files
    # --max_capus/--max_memory: Define the computing resources allowed per process
    # -profile: Use a Docker container profile
    # -work-dir: Define the working directory for Nextflow
    # -resume: Resume the pipeline if it was previously interrupted
    $NEXTFLOW_BIN run nf-core/fetchngs -r 1.12.0 \
        --input "$INPUT_CSV" \
        --outdir "$SAMPLE_DIR/fastq_data" \
        --max_memory '6.GB' \
        --max_cpus 2 \
        -profile docker \
        -work-dir "$SAMPLE_DIR/work" \
        -resume
done