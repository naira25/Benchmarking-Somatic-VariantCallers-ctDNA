#!/bin/bash

# This module fetches FASTQ files for the indicated Samples B, D and E

# Define the samples to process
samples=("B" "D" "E")

# Describe the input/output directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA/"
INPUT_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA/Samples"
NEXTFLOW_BIN="$BASE_DIR/nextflow" # Describe the path to Nextflow 

# Create a loop to obtain the FASTQ files for each sample using the nf-core/fetchngs Nextflow pipeline
for sample in "${samples[@]}"; do # Loop through each sample in the array
    echo "Processing Sample ${sample}" # Print a message indicating which sample is being processed
    
    # Create a directory for the current sample
    SAMPLE_DIR="$INPUT_DIR/Sample_${sample}" # Define the directory for the current sample
    mkdir -p "$SAMPLE_DIR"
    
    # Convert the .txt file with SRA accession ids to a .csv format required by nf-core/fetchngs
    INPUT_TXT="$INPUT_DIR/SRR_Acc_${sample}.txt" # Define the path to the input .txt file 
    INPUT_CSV="$INPUT_DIR/SRR_Acc_${sample}.csv" # Define the path to the input .csv file
    cat "$INPUT_TXT" > "$INPUT_CSV" # Convert the .txt file to .csv format (in this case, we simply copy the content as the format is compatible)

    # Run the nf-core/fetchngs Nextflow pipeline
    # Provide the input .csv file with SRA accession ids
    # Define the output directory for the downloaded FASTQ files
    # Use a container profile
    # Define the working directory for Nextflow
    # Resume the pipeline if it was previously interrupted
    $NEXTFLOW_BIN run nf-core/fetchngs \
        --input "$INPUT_CSV" \
        --outdir "$SAMPLE_DIR/fastq_data" \
        --max_memory '6.GB' \
        --max_cpus 2 \
        -profile docker \
        -work-dir "$SAMPLE_DIR/work" \
        -resume
done
