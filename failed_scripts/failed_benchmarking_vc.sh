#!/bin/bash

# ==================================================================================================================================
# FAILED BENCHMARKING SCRIPT: nf-core/variantbenchmarking PIPELINE 
# ==================================================================================================================================
# This script was developed to perform a benchmarking of variant callers using the Benchmarking 
# nf-core/variantbenchmarking pipeline from Nextflow.

# We ended up not using them due to normalization limitations. 
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"

samples=("D" "E") # Define an array with the tumor samples
subsamples=("1" "3" "5") # Define an array with the subsamples to process

# Iterate through each subsample to do the benchmarking
for sub in "${subsamples[@]}"; do
    # Iterate through each tumor sample
    for sam in "${samples[@]}"; do
    
        # Define the input CSV file required by the pipeline with information from VCF variants
        INPUT_CSV="$BASE_DIR/results/Results_Subsample_${sub}/variant_calling/benchmarking_parameters/samplesheet_${sub}_B_${sam}_benchmarking.csv"
        # Define the output directory
        OUTPUT_DIR="$BASE_DIR/results/Results_Subsample_${sub}/variant_calling/benchmarking_parameters/benchmarking_vc_B_${sam}"

        # Run the nf-core/variantbenchmarking pipeline from Nextflow
        nextflow run nf-core/variantbenchmarking -r 1.5.0 \
            -c "${BASE_DIR}/conf/nextflow.config" \
            -profile docker \
            --input "${INPUT_CSV}" \
            --outdir "${OUTPUT_DIR}" \
            --fasta "${BASE_DIR}/genome/Chr7/chr7.fa" \
            --fai "${BASE_DIR}/genome/Chr7/chr7.fa.fai" \
            --dictionary "${BASE_DIR}/genome/Chr7/chr7.dict" \
            --analysis "somatic" \
            --variant_type "small" \
            --method "sompy" \
            --truth_id "TUMOR" \
            --truth_vcf "${BASE_DIR}/metadata/truthset_files/TruthSet_chr7/KnownPositives_chr7_rename.vcf.gz" \
            --targets_bed "${BASE_DIR}/metadata/chr7_target.bed.gz" \
            --regions_bed "${BASE_DIR}/metadata/chr7_target.bed.gz" \
            --enable_missing_genotypes "test,truth" \
            --container_options "--platform linux/amd64"
    done
done