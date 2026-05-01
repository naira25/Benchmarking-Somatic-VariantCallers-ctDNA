#!/bin/bash

# ==================================================================================================================================
# TRUTH SET DATA PREPARATION 
# ==================================================================================================================================
# This script processes the Ground Truth data.

# Part1. Ground Truth variants file KnownPositives_hg19ToHg38.vcf.gz by:
# 1. Indexing the original .vcf file (.csi) to enable coordinate-based filtering
# 2. Subseting the variants to retain only those located on Chromosome 7
# 3. Standardizing the sample name by renaming it to "TUMOR" to ensure consistency with test variants files
# 4. Indexing (.csi) the processed KnownPositives_chr7_rename.vcf.gz file necessary for downstream analysis

# Part2. Ground Truth high-confidence regions file CTR_hg38.bed.gz by:
# 1. Filtering the BED file to include only High-Confidence regions within Chromosome 7
# 2. Generating an index for the processed BED file
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
TRUTH_DIR="$BASE_DIR/metadata/truthset_files"

# Create the output directory to store the processed truth set variants
OUTPUT_DIR="$TRUTH_DIR/TruthSet_chr7"
mkdir -p "$OUTPUT_DIR"

# Establish the container to pull bcftools 
BCFTOOLS_IMAGE="community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585"

# Create a file with the new sample name "TUMOR"
echo "TUMOR" > "$OUTPUT_DIR/tumor_name.txt"

# Process the KnownPositives_hg19ToHg38.vcf.gz file
docker run -v "$TRUTH_DIR":/data "$BCFTOOLS_IMAGE" sh -c "
    # 1. Index the original file
    bcftools index -f -c /data/KnownPositives_hg19ToHg38.vcf.gz && \

    # 2. Subset variants only present at chromosome 7
    bcftools view -r chr7 /data/KnownPositives_hg19ToHg38.vcf.gz -Oz -o /data/TruthSet_chr7/KnownPositives_chr7.vcf.gz && \

    # 3. Rename the truth sample to TUMOR using the previously established .txt file
    bcftools reheader -s /data/TruthSet_chr7/tumor_name.txt /data/TruthSet_chr7/KnownPositives_chr7.vcf.gz | \
    bcftools view -Oz -o /data/TruthSet_chr7/KnownPositives_chr7_rename.vcf.gz && \

    # 4. Index previously processed and renamed file
    bcftools index -c /data/TruthSet_chr7/KnownPositives_chr7_rename.vcf.gz
"

# Process the CTR_hg38.bed.gz file
docker run -v "$TRUTH_DIR":/data "$BCFTOOLS_IMAGE" sh -c "
    # Filter BED regions only present in chromosome 7
    zgrep '^chr7' /data/CTR_hg38.bed.gz | bgzip -c > /data/TruthSet_chr7/CTR_hg38_chr7.bed.gz && \
    # Generate an index
    tabix -p bed /data/TruthSet_chr7/CTR_hg38_chr7.bed.gz
"