#!/bin/bash

# ==================================================================================================================================
# VARIANT CALLERS METRICS SUMMARY
# ==================================================================================================================================
# This script generates a summary CSV file with metrics from different variant callers
# (LoFreq, VarDict, VarScan2, Mutect2 and the Ensemble Approach) for each subsample (1, 3, 5) and tumor sample (D, E) by:
# 1. Counting variants that passed all quality filters (PASS)
# 2. Categorizing identified variants into SNPs, Indels and other types
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_CSV="$BASE_DIR/results/Results_Subsample_1/variant_calling/summary.csv"

subsamples=("1" "3" "5") # Define an array with the subsamples to process
variant_callers=("lofreq" "vardict" "varscan" "mutect" "ensemble") # Define an array with the variant callers

# Create a file with the parameters to be displayed (CSV header)
echo "Subsample,Caller,Sample,PASS,SNPs,Indels,Other" > "$OUTPUT_CSV"

# Iterate through each of the subsamples variant calling results
for sub in "${subsamples[@]}"; do
    VARIANTS_DIR="$BASE_DIR/results/Results_Subsample_${sub}/variant_calling/benchmarking_vcf" # Define the path to the input .vcf files
    
    # Iterate through each of the variant callers results
    for caller in "${variant_callers[@]}"; do

        # Iterate through each of the .vcf files
        for vcf in "$VARIANTS_DIR"/*"${caller}"*.vcf.gz; do
            
            # Exclude indexed (.csi) files from the processing loop
            [[ "$vcf" == *.csi ]] && continue
            
            # Identify the tumor sample (D or E) from the filename
            filename=$(basename "$vcf")
            if [[ "$filename" == *"_D_"* ]]; then 
                mostra="D"
            else 
                mostra="E"
            fi

            # 1. Count variants with a "PASS" filter status
            pass_count=$(bcftools query -f '%FILTER\n' "$vcf" 2>/dev/null | grep -cw "PASS")
            
            # 2. Count SNPs, Indels and other variant types
            snps=$(bcftools view -v snps "$vcf" 2>/dev/null | grep -v "^#" | wc -l) # Count SNPs
            indels=$(bcftools view -v indels "$vcf" 2>/dev/null | grep -v "^#" | wc -l) # Count Indels
            altres=$(bcftools view -v other "$vcf" 2>/dev/null | grep -v "^#" | wc -l)  # Count other variants

            # Append the extracted metrics to the final CSV file
            echo "S$sub,$caller,$mostra,$pass_count,$snps,$indels,$altres" >> "$OUTPUT_CSV"
        done
    done
done