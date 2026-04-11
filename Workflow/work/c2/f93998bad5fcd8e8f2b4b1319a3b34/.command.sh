#!/bin/bash -ue
java -jar /usr/picard/picard.jar MarkDuplicates \
    -I SRX9642665_SRR13209649.aligned.sorted.bam \
    -O SRX9642665_SRR13209649_duplicates.bam \
    -M SRX9642665_SRR13209649_metrics.txt \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true

if [ -f SRX9642665_SRR13209649_duplicates.bai ]; then
    mv SRX9642665_SRR13209649_duplicates.bai SRX9642665_SRR13209649_duplicates.bam.bai
fi
