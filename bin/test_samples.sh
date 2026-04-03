#!/bin/bash

samples=("B" "D" "E")

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
TEST_DIR="$BASE_DIR/Test_Samples"

for sample in ${samples[@]}; do
    echo "Processing Sample ${sample}"

    SAMPLE_DIR="$BASE_DIR/Data/Samples/Sample_${sample}/fastq_data/fastq/" # Define the directory for the current sample
    SAMPLE_TEST_DIR="$TEST_DIR/Sample_${sample}/fastq_data/fastq/" # Define the directory for the current sample in the test directory
    
    mkdir -p "$SAMPLE_TEST_DIR" # Create a directory for the current sample in the output directory

    fastq_files=$(ls "$SAMPLE_DIR"/*.fastq.gz | sort | head -n 4)

    for fastq in $fastq_files; do
        cp "$fastq" "$SAMPLE_TEST_DIR" # Copy the current FASTQ file to the output directory
    done
done
