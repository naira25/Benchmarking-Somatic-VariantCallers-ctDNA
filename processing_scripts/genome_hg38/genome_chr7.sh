#!/bin/bash

# ==================================================================================================================================
# REFERENCE GENOME PROCESSING: Chromosome 7 Homo Sapiens (GRCh38)
# ==================================================================================================================================
# This script prepares the Chromosome 7 reference sequence from the UCSC Genome Browser.
# It performs the following steps to enable downstream alignment and variant calling by:
# 1. Downloading and decompress the raw FASTA sequence (chr7.fa.gz)
# 2. Generating BWA indexes required for the alignment process
# 3. Producing a Fasta Index (.fai) for efficient access to genomic coordinates
# 4. Creating a Sequence Dictionary (.dict) necessary for GATK and other analysis tools
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_DIR="$BASE_DIR/genome/chr7_genome"

# Establish the container to pull necessary tools
BWA_IMAGE="community.wave.seqera.io/library/bwa:0.7.19--a2905626cda4750d"
SAMTOOLS_IMAGE="community.wave.seqera.io/library/samtools:1.23.1--e8c68bc6da750dc8"

# 1. Download and decompress chromosome 7 sequence
wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
gunzip chr7.fa.gz

# 2. Build BWA index files
docker run -v $(pwd):/data "$BWA_IMAGE" \
    bwa index /data/chr7.fa

# 3. Produce FASTA index files
docker run -v $(pwd):/data "$SAMTOOLS_IMAGE" \
    samtools faidx /data/chr7.fa

# 4. Create sequence dictionary
docker run -v $(pwd):/data "$SAMTOOLS_IMAGE" \
    samtools dict /data/chr7.fa -o /data/chr7.dict