#!/bin/bash -ue
lofreq somatic \
    -n B_SRX9642665_SRR13209649_duplicates.bam \
    -t E_SRX9642659_SRR13209655_duplicates.bam \
    -f chr7.fa \
    -l chr7_target.bed \
    --threads 4 \
    -o B_SRX9642665_SRR13209649_duplicates_vs_E_SRX9642659_SRR13209655_duplicates \
    --call-indels
