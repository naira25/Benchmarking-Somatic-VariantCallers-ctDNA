#!/bin/bash

# ==================================================================================================================================
# FAILED BENCHMARKING SCRIPT: nf-core/variantbenchmarking PIPELINE 
# ==================================================================================================================================


BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
subsamples=("1" "3" "5")

variant_callers=("lofreq" "mutect" "vardict" "varscan" "ensemble")

for sub in ${subsamples[@]}; do
    
    VCF_DIR="$BASE_DIR/Results/Results_Subsample_${sub}/variant_calling/benchmarking_vcf/benchmarking_preprocessed"
    OUTPUT_DIR="$BASE_DIR/Results/Results_Subsample_${sub}/variant_calling/benchmarking_parameters"
    mkdir -p "$OUTPUT_DIR"

    CSV_FILE_D="$OUTPUT_DIR/samplesheet_${sub}_B_D_benchmarking.csv"
    CSV_FILE_E="$OUTPUT_DIR/samplesheet_${sub}_B_E_benchmarking.csv"
    
    echo "id,test_vcf,caller" > "$CSV_FILE_D"
    echo "id,test_vcf,caller" > "$CSV_FILE_E"

    for caller in "${variant_callers[@]}"; do
        
        for vcf in "$VCF_DIR"/*_D_*"$caller"*.vcf.gz; do
            [[ "$vcf" == *.csi ]] && continue

            id="TUMOR"
            echo "$id,$vcf,$caller" >> "$CSV_FILE_D"
        done

        for vcf in "$VCF_DIR"/*_E_*"$caller"*.vcf.gz; do
            [ -e "$vcf" ] || continue
            [[ "$vcf" == *.csi ]] && continue

            id="TUMOR"
            echo "$id,$vcf,$caller" >> "$CSV_FILE_E"
        done
        
    done
done