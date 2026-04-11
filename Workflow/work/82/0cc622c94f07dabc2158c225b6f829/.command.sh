#!/bin/bash -ue
samtools sort -o "SRX9642659_SRR13209655.aligned.sorted.bam" SRX9642659_SRR13209655.aligned.raw.bam
samtools index "SRX9642659_SRR13209655.aligned.sorted.bam"
