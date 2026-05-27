#!/bin/bash

# ==================================================================================================================================
# PREPROCESSING VCF FILES FOR BENCHMARKING
# ==================================================================================================================================
# This script standardizes VCF files from different callers (Mutect2, VarDict, VarScan2, LoFreq and the Ensemble approach) by:
# 1. Extracting only the tumor sample from multi-sample VCFs.
# 2. Renaming the tumor sample column to a generic "TUMOR" header.
# 3. For LoFreq, it adds missing genotype (GT) information to ensure compatibility.
# 4. Indexing all resulting files for downstream analysis.
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"

# Define the subsamples and their specific SRA accessions (Subsample:Normal:TumorD:TumorE)
subsamples_str=(
    "1:SRX9642665_SRR13209649:SRX9642709_SRR13209605:SRX9642659_SRR13209655"
    "3:SRX9642673_SRR13209641:SRX9642677_SRR13209637:SRX9642685_SRR13209629"
    "5:SRX9642700_SRR13209614:SRX9642703_SRR13209611:SRX9642712_SRR13209602"
)

# Create a temporary file with the new sample name "TUMOR"
echo "TUMOR" > tumor_header.txt

# Iterate through each subsample files
for structure in "${subsamples_str[@]}"; do

    # Split the strings into variables
    IFS=":" read -r subsample normal tumorD tumorE <<< "$structure"
    
    # Define input VCF files for each variant caller
    INPUT_DIR="${BASE_DIR}/results/Results_Subsample_${subsample}/variant_calling/benchmarking_vcf"
    # Define output processed VCF files
    OUTPUT_DIR="${BASE_DIR}/results/Results_Subsample_${subsample}/variant_calling/benchmarking_vcf/benchmarking_preprocessed"
    mkdir -p "$OUTPUT_DIR"

    samples=("D" "E") # Define an array with the tumor samples to process
    # Iterate through each tumor sample
    for sam in ${samples[@]}; do

    # Assign the correct SRR accession based on the current sample (D or E)
    [[ "$sam" == "D" ]] && tumor_sample="$tumorD" || tumor_sample="$tumorE"
    normal_sample="$normal"
    
        # Process Mutect2 VCFs
        # Extract tumor sample and rename column to "TUMOR"
        bcftools view -s "$tumor_sample" "${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_mutect_final.vcf.gz" | \
        bcftools reheader -s tumor_header.txt | \
        bcftools view -Oz -o "${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_mutect_tumor.vcf.gz"

        # Process VarDict VCFs
        # Extract tumor sample and rename column to "TUMOR"
        bcftools view -s "${sam}_${tumor_sample}" "${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_vardict_final.vcf.gz" | \
        bcftools reheader -s tumor_header.txt | \
        bcftools view -Oz -o "${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_vardict_tumor.vcf.gz"

        # Process VarScan2 VCFs
        # Extract already existing TUMOR column
        bcftools view -s "TUMOR" ${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_varscan_final.vcf.gz -Oz -o ${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_varscan_tumor.vcf.gz

        # Process LoFreq VCFs
        # Add genotype (GT) column
        gzcat ${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_lofreq_final.vcf.gz | grep "^##" > header_nou.vcf
        echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> header_nou.vcf
        # Append GT column/values to data lines and combine with new header
        gzcat ${INPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_lofreq_final.vcf.gz | grep -v "^##" | \
        awk 'BEGIN {FS="\t"; OFS="\t"} /^#CHROM/ {print $0, "FORMAT", "TUMOR"; next} {print $0, "GT", "0/1"}' | \
        cat header_nou.vcf - | \
        bcftools view -Oz -o ${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_lofreq_tumor.vcf.gz

        # Process Ensemble approach VCFs
        bcftools view -s "TUMOR" ${INPUT_DIR}/S${subsample}_B_${sam}_ensemble_2.vcf.gz -Oz -o ${OUTPUT_DIR}/S${subsample}_B_${normal}_vs_${sam}_${tumor_sample}_ensemble_tumor.vcf.gz
    done

    # Index all processed tumor-only VCF files
    for vcf in ${OUTPUT_DIR}/*.vcf.gz; do
        [[ "$vcf" == *.csi ]] && continue
        bcftools index -c "$vcf" 
    done
done