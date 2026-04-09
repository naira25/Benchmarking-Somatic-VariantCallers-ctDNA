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
params.reads = "${params.project}/Data/Sample_*/fastq_data/fastq/*_{1,2}.fastq.gz"
params.genome ="${params.project}/Data/genome"
params.outdir = "${params.project}/Data/Results"


/*
 * Workflow Channels
 */

reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
genome_ch = Channel.fromPath(params.genome, checkIfExists: true)

/*
 * Workflow Processes
 */

process FASTQC {

    tag{"QC ${sample_id}"}

    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fasqc_results", emit: fastqc_reports

    script:
    """
    mkdir fasqc_results
    fastqc -o fasqc_results ${reads} -t 2 -q
    """ 
}

process BWAMEM_SAMTOOLS_SORT {

    tag{"Alinging and Sorting ${sample_id}"}

    container 'biocontainers/bwa:v0.7.17-3-deb_cv1'
 
    input:
    tuple val(sample_id), path(reads)
    path(indexed_reference_genome)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), emit: bam_alingments_sorted
    
    script:
    """
    bwa mem -t -2 -R ${indexed_reference_genome} ${reads[0]} ${reads[1]} | \\
    samtools view -bS - | \\
    samtools sort -o "${sample_id}.aligned.sorted.bam" -
    """
}

process PICARD_MARKDUPLICATES {

    tag{"Marking Duplicates_${sample_id}"}

    container 'biocontainers/picard:v1.139_cv3'

    input:
    tuple val(sample_id), path(sorted_alingment)

    output:
    tuple val(sample_id), path("${sample_id}_duplicates.bam"), emit: duplicates_alignments
    path "${sample_id}_metrics.txt", emit: metrics

    script:
    """
    java -jar picard.jar MarkDuplicates \\
        I=${sorted_alingment} \\
        O=${sample_id}_duplicates.bam \\
        M=${sample_id}_metrics.txt \\
        REMOVE_DUPLICATES=false
    """ 
}


/*
 * Workflow
 */

workflow {

    main:
        FASTQC(reads_ch)
        BWAMEM_SAMTOOLS_SORT(reads_ch, genome_ch)
        PICARD_MARKDUPLICATES(BWAMEM_SAMTOOLS_SORT.out)

    publish:
        fastqc_out = FASTQC.out
        alingments_sorted_out = BWAMEM_SAMTOOLS_SORT.out
        duplicates_out = PICARD_MARKDUPLICATES.out.duplicates_alignments
        metrics_out = PICARD_MARKDUPLICATES.out.metrics
}

output {
    fastqc_out {
        path "${params.outdir}/fastqc_reports"
        mode 'copy'
    }
    alingments_sorted {
        path "${params.outdir}/bam_alingments_sorted"
        mode 'copy'
    }
    duplicates_out {
        path "${params.outdir}/mark_duplicates/duplicates_alignments"
        mode 'copy'
    }
    metrics_out {
        path "${params.outdir}/mark_duplicates/metrics"
        mode 'copy'
    }
}