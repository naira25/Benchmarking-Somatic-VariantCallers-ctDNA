
samples=("B" "D" "E")

INPUT_DIR="/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA/SamplesSample_*/fastq_data/fastq"

for sample in ${samples[@]}; do
    
    Patient=
    Sample=
    Lane="Lane"
    fastq_1=INPUT_DIR/Sample_*/fastq_data/fastq/*_1.fastq.gz
    fastq_2=INPUT_DIR/Sample_*/fastq_data/fastq/*_1.fastq.gz