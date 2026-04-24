#!/bin/bash

BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_CSV="$BASE_DIR/Results/Results_Subsample_1/variant_calling/summary.csv"

subsamples=("1" "3" "5")
variant_callers=("lofreq" "vardict" "varscan" "mutect" "ensemble")

echo "Subsample,Caller,Sample,PASS,SNPs,Indels,Other" > "$OUTPUT_CSV"

for sub in "${subsamples[@]}"; do
    VARIANTS_DIR="$BASE_DIR/Results/Results_Subsample_${sub}/variant_calling/benchmarking_vcf"
    
    for caller in "${variant_callers[@]}"; do

        for vcf in "$VARIANTS_DIR"/*"${caller}"*.vcf.gz; do
            
            [[ "$vcf" == *.csi ]] && continue
            
            filename=$(basename "$vcf")
            
            if [[ "$filename" == *"_D_"* ]]; then 
                mostra="D"
            else 
                mostra="E"
            fi

            pass_count=$(bcftools query -f '%FILTER\n' "$vcf" 2>/dev/null | grep -cw "PASS")
            snps=$(bcftools view -v snps "$vcf" 2>/dev/null | grep -v "^#" | wc -l)
            indels=$(bcftools view -v indels "$vcf" 2>/dev/null | grep -v "^#" | wc -l)
            altres=$(bcftools view -v other "$vcf" 2>/dev/null | grep -v "^#" | wc -l)

            echo "S$sub,$caller,$mostra,$pass_count,$snps,$indels,$altres" >> "$OUTPUT_CSV"
        done
    done
done