#!/bin/bash

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_DIR="$BASE_DIR/Genome/Chr7"

echo "Downloading chr7 reference genome from UCSC"
wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
gunzip chr7.fa.gz

echo "Indexing reference genome with BWA"
docker run --rm -v $(pwd):/data community.wave.seqera.io/library/bwa:0.7.19--a2905626cda4750d \
    bwa index /data/chr7.fa

echo "Indexing reference genome with Samtools"
docker run --rm -v $(pwd):/data community.wave.seqera.io/library/samtools:1.23.1--e8c68bc6da750dc8 \
    samtools faidx /data/chr7.fa

echo "Creating .dict dictionary"
docker run --rm -v $(pwd):/data community.wave.seqera.io/library/samtools:1.23.1--e8c68bc6da750dc8 \
    samtools dict /data/chr7.fa -o /data/chr7.dict