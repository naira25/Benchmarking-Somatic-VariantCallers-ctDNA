#!/bin/bash

BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
ALINGMENTS_DIR="$BASE_DIR/Results/mark_duplicates/duplicates_alignments"
TXT_FILE="${BASE_DIR}/Samples"

for alignment in "$ALINGMENTS_DIR"/*.{bam,bai}; do
    
    file_name=$(basename "$alignment")

    srr=$(echo "$file_name" | cut -d'_' -f2)

    sample=""
    
    if grep -q "$srr" "$TXT_FILE/SRR_Acc_B.txt"; then
        sample="B"
    
    elif grep -q "$srr" "$TXT_FILE/SRR_Acc_D.txt"; then
        sample="D"
    
    elif grep -q "$srr" "$TXT_FILE/SRR_Acc_E.txt"; then
        sample="E"
    fi

    if [ -n "$sample" ]; then
        mv "$alignment" "$ALINGMENTS_DIR/${sample}_${file_name}"
    fi
done