#!/bin/bash

samples=("B" "D" "E")
subsamples=("1" "2" "3" "4" "5" "6")
labs=("10" "20" "21")

BASE_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
OUTPUT_DIR="$BASE_DIR/Subsamples"

# Crear directoris
for sub in "${subsamples[@]}"; do
    mkdir -p "$OUTPUT_DIR/Subsample_${sub}"
done

for sample in "${samples[@]}"; do

    SAMPLE_DIR="$BASE_DIR/Data/Samples/Sample_${sample}/fastq_data/fastq"
    CSV_FILE="${SAMPLE_DIR}/../samplesheet/samplesheet.csv"
    TXT_FILE="${SAMPLE_DIR}/../samplesheet/samplesheet.txt"

    # Neteja de format: canviem només les comes que separen columnes (" , ") per tabuladors
    # Això garanteix que la columna 22 segueixi sent la de les bases
    sed 's/","/\t/g; s/"//g' "$CSV_FILE" > "$TXT_FILE"

    for lab in "${labs[@]}"; do
        if [ "$lab" == "10" ]; then base_folder=1; fi
        if [ "$lab" == "20" ]; then base_folder=3; fi
        if [ "$lab" == "21" ]; then base_folder=5; fi

        # Obtenim els IDs de les 2 millors rèpliques
        top_ids=$(grep "ST$lab" "$TXT_FILE" | sort -t$'\t' -k22,22rn | head -n 2 | cut -f4)

        i=0 # Reiniciem el comptador per a cada laboratori
        for id in $top_ids; do
            target_sub=$((base_folder + i))
            cp "$SAMPLE_DIR"/*"${id}"*.fastq.gz "$OUTPUT_DIR/Subsample_${target_sub}/"
            ((i++)) 
        done
    done
    rm "$TXT_FILE"
done