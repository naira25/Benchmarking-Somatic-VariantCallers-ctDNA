/*
==============================================
    Nextflow Workflow for Variant Calling
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

process MUTECT2 {

    tag{"Variant Calling ${sample_id} with Mutect2"}

    container 'broadinstitute/gatk'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fasqc_results", emit: fastqc_reports

    script:
    """
    gatk Mutect2 \
        -I tumor.bam \
        -I normal.bam \
        -O somatic.vcf.gz \
        -R reference.fa \
        --normal-sample normal_sample_name

    gatk FilterMutectCalls \
        -R reference.fa \
        -V ${sample_id}_unfiltered.vcf.gz \
        -O ${sample_id}_filtered.vcf.gz
    """ 
}

process LOFREQ {

    tag{"Variant Calling ${sample_id} with LoFreq"}

    container 'broadinstitute/gatk'
 
    input:
    tuple val(sample_id), path(reads)
    path(indexed_reference_genome)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), emit: bam_alingments_sorted
    
    script:
    """
    lofreq somatic \
        -n normal.bam -t tumor.bam -f hg19.fa \
        --threads 8 -o out_ [-d dbsnp.vcf.gz] \
        --call-indels
    """
}

process VARDICT {

    tag{"Variant Calling ${sample_id} with Vardict"}

    container 'broadinstitute/gatk'

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

    process ABEMUS {

    tag{"Variant Calling ${sample_id} with ABEMUS"}

    container 'broadinstitute/gatk'

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