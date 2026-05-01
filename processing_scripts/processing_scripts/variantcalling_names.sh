#!/bin/bash

# ==================================================================================================================================
# ALIGNMENTS RENAME
# ==================================================================================================================================
# This script renames sorted alignments with Marked Duplicated (.bam .bai) by adding the corresponding sample B, D or E on the filename
# It identifies the sample type by cross-referencing the SRR accession number with the provided metadata text files.
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
TXT_FILE="${BASE_DIR}/samples"

subsamples=("3" "5") # Define an array with the subsamples to process

# Iterate through each of the samples alignment results
for sub in "${subsamples[@]}"; do
    ALIGNMENTS_DIR="$BASE_DIR/results/Results_Subsample_${sub}/mark_duplicates/duplicates_alignments" # Define the path to the input alignment files

    # Iterate through each alignment file (.bam and .bai) in the directory
    for alignment in "$ALIGNMENTS_DIR"/*.{bam,bai}; do
    
        # Extract the base filename for each file
        file_name=$(basename "$alignment")
        # Extract the SRR ID for each filename
        srr=$(echo "$file_name" | cut -d'_' -f2)

        # Initialize the sample variable
        sample=""
        #  Identify the SRR IDs on the accession lists to identify the sample type
        if grep -q "$srr" "$TXT_FILE/SRR_Acc_B.txt"; then
            sample="B"
        
        elif grep -q "$srr" "$TXT_FILE/SRR_Acc_D.txt"; then
            sample="D"
        
        elif grep -q "$srr" "$TXT_FILE/SRR_Acc_E.txt"; then
            sample="E"
        fi

        # Add sample name in the alignments according to theire present in a .txt file
        if [ -n "$sample" ]; then
            mv "$alignment" "$ALIGNMENTS_DIR/${sample}_${file_name}"
        fi
    done
done