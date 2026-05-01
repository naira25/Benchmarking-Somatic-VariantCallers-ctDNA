#!/bin/bash

# ==================================================================================================================================
# VARIANT CALLING ENSEMBLE APPROACH
# ==================================================================================================================================
# This script integrates results from LoFreq, Mutect2, VarScan, and VarDict. For each sample pair (B-D and B-E)
# it implements a consensus approach to retain variants identified by at least 2 out of the 4 variant callers by:
# 1. Identifying and extract variants supported by a minimum of 2 out of the 4 callers.
# 2. Merging consensus variants into a single Ensemble .vcf file and sort by genomic coordinate order 
# 3. Generating indexes (.csi) necessary for downstream benchmarking analysis
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"

subsamples=("1" "3" "5") # Define an array with the subsamples to process
samples=("D" "E") # Define an array with the tumor samples to process

# Iterate through each of the subsamples variant calling results
for sub in "${subsamples[@]}"; do

    VARIANTS_DIR="$BASE_DIR/results/Results_Subsample_${sub}/variant_calling/benchmarking_vcf" # Define the path to the input .vcf files
    cd "$VARIANTS_DIR" || continue # Move to the input .vcf files directory

    # Iterate through each of the sample pairs B-D or B-E to create the ensemble
    for sam in "${samples[@]}"; do

        # 1. Extract variants present at least in 2 variant callers and put intersect results in a new directory
        bcftools isec -n+2 -c all *_${sam}_*final.vcf.gz -p "ensemble_B_${sam}_2"
        # Iterate through each intersection file to index each compressed file
        for intersection in ensemble_B_${sam}_2/000*.vcf; do
            bgzip -f "$intersection"
            bcftools index -f "${intersection}.gz"
        done

        # 2. Merge intersection results and sort by genomic coordinate order
        bcftools merge --force-samples -m none ensemble_B_${sam}_2/000*.vcf.gz | bcftools sort -Oz -o "S${sub}_B_${sam}_ensemble_2.vcf.gz"

        # 3. Index the final Ensemble VCF file
        bcftools index -c "S${sub}_B_${sam}_ensemble_2.vcf.gz"
        
    done

    # Return to the base directory for the next subsample iteration
    cd "$BASE_DIR"
done