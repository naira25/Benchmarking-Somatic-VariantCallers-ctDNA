#!/bin/bash

# ==================================================================================================================================
# SUBSAMPLES CREATION
# ==================================================================================================================================
# This script creates 3 subsamples (1, 3 and 5) containing FASTQ files for samples B, D and E.
# It selects runs with the highest number of sequenced bases from three sequencing labs:
# - Lab 10 (corresponding to Subsample 1).
# - Lab 20 (corresponding to Subsample 3).
# - Lab 21 (corresponding to Subsample 5).
# ==================================================================================================================================

# Describe the working directories
BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_DIR="$BASE_DIR/subsamples"

samples=("B" "D" "E") # Define an array with the samples to process
subsamples=("1" "3" "5") # Define an array with the subsamples to process
labs=("10" "20" "21") # Define an array with the sequencing labs

# Iterate to create the output directory to store each subsamples FASTQ files
for sub in "${subsamples[@]}"; do
    mkdir -p "$OUTPUT_DIR/subsample_${sub}"
done

# Iterate to copy corresponding FASTQ files for samples B, D and E with a higher number of sequenced bases from each sequencing lab 10, 20 and 21 in the Subsamples directories
for sample in "${samples[@]}"; do

    # Describe the working directories
    SAMPLE_DIR="$BASE_DIR/samples/Sample_${sample}/fastq_data/fastq" # Define the path to the input FASTQ file
    CSV_FILE="${SAMPLE_DIR}/../samplesheet/samplesheet.csv" # Define the path to the input .csv file with sample information (SRA ID, lab and number of sequenced bases)
    TXT_FILE="${SAMPLE_DIR}/../samplesheet/samplesheet.txt" # Define the path to the input .txt file with sample information (SRA ID, lab and number of sequenced bases)

    # Convert the CSV samplesheet with information from each sample to a TXT file
    sed 's/","/\t/g; s/"//g' "$CSV_FILE" > "$TXT_FILE"

    # Iterate through samples from each lab to obtain the ones with a higher number of sequenced bases
    for lab in "${labs[@]}"; do
        if [ "$lab" == "10" ]; then base_folder=1; fi # Define Subsample 1, corresponding to FASTQ files for samples B, D and E with a higher number of sequenced bases in sequencing lab 10
        if [ "$lab" == "20" ]; then base_folder=3; fi # Define Subsample 3, corresponding to FASTQ files for samples B, D and E with a higher number of sequenced bases in sequencing lab 20
        if [ "$lab" == "21" ]; then base_folder=5; fi # Define Subsample 5, corresponding to FASTQ files for samples B, D and E with a higher number of sequenced bases in sequencing lab 21

        # Obtain SRA accession IDs for samples from each lab with a higher number of sequenced bases
        higher_seq_file=$(grep "ST$lab" "$TXT_FILE" | sort -t$'\t' -k22,22rn | head -n 2 | cut -f4)

        # Restart the process for each sample
        i=0
        # Iterate through previously filtered samples
        for file in $higher_seq_file; do
            target_sub=$((base_folder + i)) # Calculate the destination folder index based on the lab's base folder
            cp "$SAMPLE_DIR"/*"${file}"*.fastq.gz "$OUTPUT_DIR/Subsample_${target_sub}/" # Copy the FASTQC files from filtered samples in the Subsample directories
            ((i++)) 
        done
    done
    rm "$TXT_FILE" # Remove the temporarily created .txt file
done