#!/bin/bash -ue
lofreq somatic \
    -n B_SRX9642665_SRR13209649_duplicates.bam \
    -t D_SRX9642709_SRR13209605_duplicates.bam \
    -f chr7.fa \
    -l chr7_target.bed \
    --threads 4 \
    -o B_SRX9642665_SRR13209649_duplicates_vs_D_SRX9642709_SRR13209605_duplicates \
    --call-indels
