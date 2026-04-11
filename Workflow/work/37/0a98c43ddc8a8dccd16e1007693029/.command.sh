#!/bin/bash -ue
samtools sort -o "SRX9642665_SRR13209649.aligned.sorted.bam" SRX9642665_SRR13209649.aligned.raw.bam
samtools index "SRX9642665_SRR13209649.aligned.sorted.bam"
