/*
==================================================================================================================================
# NEXTFLOW PREPROCESSING PIPELINE
==================================================================================================================================
# This script processes paired-end raw FASTQ reads from samples B, D, and E following GATK4 best practices:
# 1. Quality Control (FastQC): generate quality reports for raw FASTQ reads using FastQC
# 2. Genome Alignment (.bam/.bai): align reads to the Chromosome 7 reference (GRCh38) sequence using BWA-MEM
#    Raw alignments are not stored locally to optimize space
# 3. Post-alignment: sort and index the resulting BAM files using Samtools
# 4. Duplicate Marking: identify PCR duplicates using Picard MarkDuplicates
#
# For each step, a specific Seqera container (compatible with bioconda and arm64x) has been used.
#
# To manage hardware limitations, the pipeline is designed to be implemented separately for each sample group
# to ensure stability and efficient resource consumption on local machines.
# ==================================================================================================================================
*/

nextflow.enable.dsl=2

/*
 * Workflow Input Parameters
 */
params.project = "/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
// Path to raw paired-end {1,2} FASTQ reads
params.reads = "${params.project}/Subsamples/Subsample_1/*_{1,2}.fastq.gz" 
// Path to Chromosome 7 reference FASTA and its corresponding index files
params.genome_sequence = "${params.project}/Genome/Chr7/chr7.fa" 
// Path to the main Results directory
params.outdir = "${params.project}/Results"


/*
 * Workflow Processes
 */

// 1. Quality Control Process using FastQC
process fastQC {

    tag "Quality Control for ${sample_id} sample"
    container 'community.wave.seqera.io/library/fastqc:0.12.1--df99cb252670875a'
    publishDir "${params.outdir}/fastqc_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc"

    script:
    """
    # Create output directory for the specific sample
    mkdir ${sample_id}_fastqc
    # Run FastQC: -o (output dir), -t (threads), -q (quiet mode)
    fastqc -o ${sample_id}_fastqc ${reads} -t 2 -q
    """ 
}

// 2. Genome Alignment to Chr7 using BWA-MEM
process bwamemAlignment {

    tag "Aligning sample ${sample_id}"
    container 'community.wave.seqera.io/library/bwa:0.7.19--a2905626cda4750d'
    
    input:
    tuple val(sample_id), path(reads)
    path indexed_files 

    output:
    tuple val(sample_id), path("${sample_id}.aligned.raw.bam")
    
    script:
    // Identify the genome .fa file among the provided index files
    def fasta = indexed_files.find { it.name.endsWith(".fa") }
    """
    # Map reads using BWA-MEM with Read Group information (-R) and redirect to BAM
    bwa mem -t 2 -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
    ${fasta} ${reads[0]} ${reads[1]} > ${sample_id}.aligned.raw.bam
    """
}

// 3. Sort and Index alignments using Samtools
process samtoolsSort {

    tag "Sorting sample ${sample_id} alignments"
    container 'community.wave.seqera.io/library/samtools:1.23.1--e8c68bc6da750dc8'
    publishDir "${params.outdir}/bam_alignments_sorted", mode: 'copy'
 
    input:
    tuple val(sample_id), path(raw_bam)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), path("${sample_id}.aligned.sorted.bam.bai")
    
    script:
    """
    # Sort the BAM file by coordinate order
    samtools sort -o "${sample_id}.aligned.sorted.bam" ${raw_bam}
    # Create a coordinate index (.bai) for the sorted BAM
    samtools index "${sample_id}.aligned.sorted.bam"
    """
}

// 4. Mark PCR duplicates on sorted alignments using Picard MarkDuplicates
process picardMarkDuplicates {

    tag "Marking Duplicates ${sample_id}"
    container 'community.wave.seqera.io/library/picard:3.4.0--6f28fdc142d7e8d3'
    publishDir "${params.outdir}/mark_duplicates/duplicates_alignments", mode: 'copy', pattern: "*.bam*"
    publishDir "${params.outdir}/mark_duplicates/metrics", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample_id), path(sorted_alignment)

    output:
    tuple val(sample_id), path("${sample_id}_duplicates.bam"), path("${sample_id}_duplicates.bam.bai"), emit: duplicates_alignments
    path "${sample_id}_metrics.txt", emit: metrics

    script:
    """
    # Identify PCR duplicates and generate a metrics file
    picard MarkDuplicates \\
        -I ${sorted_alignment} \\
        -O ${sample_id}_duplicates.bam \\
        -M ${sample_id}_metrics.txt \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true

    # Standardize the index file name to .bam.bai if Picard created it as .bai
    if [ -f ${sample_id}_duplicates.bai ]; then
        mv ${sample_id}_duplicates.bai ${sample_id}_duplicates.bam.bai
    fi
    """ 
}


/*
 * Main Workflow
 */
workflow {

    // Channel for paired-end raw FASTQ reads
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    
    // Channel for Chromosome 7 reference FASTA
    genome_ch = Channel.fromPath("${params.genome_sequence}*", checkIfExists: true).collect()

    // Processes execution
    fastQC(reads_ch) // 1. Quality Control
    bwa_raw_ch = bwamemAlignment(reads_ch, genome_ch) // 2. Alignment
    bam_sorted_ch = samtoolsSort(bwa_raw_ch) // 3. Sorting and Indexing Alignments
    picardMarkDuplicates(bam_sorted_ch.map{ it -> [it[0], it[1]] }) // 4. Mark Duplicates
}