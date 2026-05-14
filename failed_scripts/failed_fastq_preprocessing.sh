#!/bin/bash

# ==================================================================================================================================
# FAILED PREPROCESSING SCRIPT: nf-core/sarek PIPELINE 
# ==================================================================================================================================
# This script was originally developed to preprocess raw FASTQ reads for samples B, D, and E 
# using the nf-core/sarek Nextflow pipeline.

# At first, we tested the pipeline using the whole Homo Sapiens reference genome (GRCh38). As the process was computationally
# consuming, we narrowed the analysis and focused only on Chromosome 7. However, the pipeline's resource demands still exceeded the 
# capabilities of the local infrastructure. Therefore, we developed our own preprocessing workflow using nf-core/sarek as a reference.
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"

subsamples=("1" "3" "5") # Define an array with the subsamples to process

# Iterate through each subsample to do the benchmarking
for sub in "${subsamples[@]}"; do

    # Define the input CSV file required by the pipeline with information from FASTQ files
    INPUT_CSV="$BASE_DIR/metadata/samplesheet.csv"
    # Define the output directory
    OUTPUT_DIR="$BASE_DIR/results/Results_Subsample_*"

    # Run the nf-core/sarek pipeline from Nextflow
    nextflow run nf-core/sarek -r 3.8.1 \
        -c ${BASE_DIR}/conf/nextflow.config \
        -profile docker \
        -resume \
        --input ${BASE_DIR}/metadata/samplesheet.csv \
        --outdir ${OUTPUT_DIR} \
        --genome GATK.GRCh38 \
        --wes true \
        --intervals ${BASE_DIR}/metadata/ROC2.bed \
        --save_output_as_bam true \
        --save_reference false \
        --skip_tools baserecalibrator
done