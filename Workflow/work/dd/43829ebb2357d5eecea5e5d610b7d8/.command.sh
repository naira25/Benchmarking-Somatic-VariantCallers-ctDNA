#!/bin/bash -ue
samtools sort -o "SRX9642709_SRR13209605.aligned.sorted.bam" SRX9642709_SRR13209605.aligned.raw.bam
samtools index "SRX9642709_SRR13209605.aligned.sorted.bam"
