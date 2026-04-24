/*
==============================================
    Nextflow Workflow for fastq quality control, alignment, sorting and indexing, marking duplicates
===============================================
*/

nextflow.enable.dsl=2

/*
 * Workflow Input Parameters
 */

params.project = "/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA"
params.reads = "${params.project}/Subsamples/Subsample_1/*_{1,2}.fastq.gz"
params.genome_sequence = "${params.project}/Genome/references/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
params.outdir = "${params.project}/Results"


/*
 * Workflow Channels
 */

reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
genome_ch = Channel.fromPath("${params.genome_sequence}*", checkIfExists: true).collect()

/*
 * Workflow Processes
 */

process FASTQC {

    tag{"QC ${sample_id}"}

    container 'biocontainers/fastqc:v0.11.9_cv8'

    publishDir "${params.outdir}/fastqc_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc"

    script:
    """
    mkdir ${sample_id}_fastqc
    fastqc -o ${sample_id}_fastqc ${reads} -t 2 -q
    """ 
}

process BWAMEM_ALIGNMENT {

    tag{"Alinging ${sample_id}"}

    container 'biocontainers/bwa:v0.7.17-3-deb_cv1'
 
    input:
    tuple val(sample_id), path(reads)
    path indexed_files 

    output:
    tuple val(sample_id), path("${sample_id}.aligned.raw.bam")
    
    script:
    def fasta = indexed_files.find { it.name.endsWith(".fa") }
    """
    bwa mem -t 2 -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
    ${fasta} ${reads[0]} ${reads[1]} > ${sample_id}.aligned.raw.bam
    """
}

process SAMTOOLS_SORT {

    tag{"Sorting ${sample_id}"}

    container 'biocontainers/samtools:v1.9-4-deb_cv1'

    publishDir "${params.outdir}/bam_alingments_sorted", mode: 'copy'
 
    input:
    tuple val(sample_id), path(raw_bam)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), path("${sample_id}.aligned.sorted.bam.bai")
    
    script:
    """
    samtools sort -o "${sample_id}.aligned.sorted.bam" ${raw_bam}
    samtools index "${sample_id}.aligned.sorted.bam"
    """
}

process PICARD_MARKDUPLICATES {

    tag{"Marking Duplicates_${sample_id}"}

    container 'biocontainers/picard:v1.139_cv3'

    publishDir "${params.outdir}/mark_duplicates/duplicates_alignments", mode: 'copy', pattern: "*.bam*"
    publishDir "${params.outdir}/mark_duplicates/metrics", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample_id), path(sorted_alingment)

    output:
    tuple val(sample_id), path("${sample_id}_duplicates.bam"), path("${sample_id}_duplicates.bam.bai"), emit: duplicates_alignments
    path "${sample_id}_metrics.txt", emit: metrics

    script:
    """
    java -jar /usr/local/bin/picard.jar MarkDuplicates \\
        I=${sorted_alingment} \\
        O=${sample_id}_duplicates.bam \\
        M=${sample_id}_metrics.txt \\
        REMOVE_DUPLICATES=false
    
    samtools index ${sample_id}_duplicates.bam
    """ 
}


/*
 * Workflow
 */

workflow {

        FASTQC(reads_ch)
        bwa_raw_ch = BWAMEM_ALIGNMENT(reads_ch, genome_ch)
        bam_sorted_ch = SAMTOOLS_SORT(bwa_raw_ch)
        PICARD_MARKDUPLICATES(bam_sorted_ch.map{ it -> [it[0], it[1]] })
}