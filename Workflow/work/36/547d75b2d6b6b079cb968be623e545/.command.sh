#!/bin/bash -ue
java -jar /usr/picard/picard.jar MarkDuplicates \
    -I SRX9642709_SRR13209605.aligned.sorted.bam \
    -O SRX9642709_SRR13209605_duplicates.bam \
    -M SRX9642709_SRR13209605_metrics.txt \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true

if [ -f SRX9642709_SRR13209605_duplicates.bai ]; then
    mv SRX9642709_SRR13209605_duplicates.bai SRX9642709_SRR13209605_duplicates.bam.bai
fi
