#!/bin/bash


# ==================================================================================================================================
# PREPROCESSING VCF FILES FOR BENCHMARKING
# ==================================================================================================================================
# This script standardizes VCF files from different callers (MuTect2, VarDict, VarScan, LoFreq, Ensemble) by:
# 1. Extracts only the tumor sample from multi-sample VCFs.
# 2. Renames the tumor sample column to a generic "TUMOR" header.
# 3. For LoFreq, it adds missing genotype (GT) information to ensure compatibility.
# 4. Indexes all resulting files for downstream analysis.
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"

# Define the subsamples and their specific SRR accessions (Subsample:Normal:TumorD:TumorE)
subsamples_str=(
    "1:SRX9642665_SRR13209649:SRX9642709_SRR13209605:SRX9642659_SRR13209655"
    "3:SRX9642673_SRR13209641:SRX9642677_SRR13209637:SRX9642685_SRR13209629"
    "5:SRX9642700_SRR13209614:SRX9642703_SRR13209611:SRX9642712_SRR13209602"
)

# Create a temporary file with the new sample name "TUMOR"
echo "TUMOR" > tumor_header.txt

for structure in "${subsamples_str[@]}"; do

    IFS=":" read -r subsample normal tumorD tumorE <<< "$structure"

    INPUT_DIR="${BASE_DIR}/results/Results_Subsample_${subsample}/variant_calling/benchmarking_vcf"
    OUTPUT_DIR="${BASE_DIR}/results/Results_Subsample_${subsample}/variant_calling/benchmarking_vcf/benchmarking_preprocessed"
    mkdir -p "$OUTPUT_DIR"

    samples=("D" "E")
    for sam in ${samples[@]}; do

    [[ "$sam" == "D" ]] && tumor_sample="$tumorD" || tumor_sample="$tumorE"
    normal_sample="$normal"
    
        bcftools view -s "$tumor_sample" "${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_mutect_final.vcf.gz" | \
        bcftools reheader -s tumor_header.txt | \
        bcftools view -Oz -o "${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_mutect_tumor.vcf.gz"

        # --- VARDICT (ID intern sol portar D_) ---
        bcftools view -s "${sam}_${tumor_sample}" "${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_vardict_final.vcf.gz" | \
        bcftools reheader -s tumor_header.txt | \
        bcftools view -Oz -o "${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_vardict_tumor.vcf.gz"

        # --- VARSCAN ---
        bcftools view -s "TUMOR" ${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_varscan_final.vcf.gz -Oz -o ${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_varscan_tumor.vcf.gz

        # --- LOFREQ ---
        gzcat ${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_lofreq_final.vcf.gz | grep "^##" > header_nou.vcf
        echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> header_nou.vcf

        # 2. Ajuntem el header nou amb les dades modificades per awk i comprimim
        gzcat ${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_lofreq_final.vcf.gz | grep -v "^##" | \
        awk 'BEGIN {FS="\t"; OFS="\t"} /^#CHROM/ {print $0, "FORMAT", "TUMOR"; next} {print $0, "GT", "0/1"}' | \
        cat header_nou.vcf - | \
        bcftools view -Oz -o ${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_lofreq_tumor.vcf.gz

        # 3. Netegem
        rm header_nou.vcf

        # --- ENSEMBLE ---
        # Provem amb el SRR i si falla amb "TUMOR"
        bcftools view -s "TUMOR" ${INPUT_DIR}/S${subsample}_B_${sam}_ensemble_2.vcf.gz -Oz -o ${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_ensemble_tumor.vcf.gz
    done

    for vcf in ${OUTPUT_DIR}/*.vcf.gz; do
        [[ "$vcf" == *.csi ]] && continue
        bcftools index -c "$vcf" 
    done
done

# Remove the temporarily created file
rm tumor_header.txt