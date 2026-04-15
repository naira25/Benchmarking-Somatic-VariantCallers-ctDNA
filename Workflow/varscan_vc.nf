/*
==============================================
    Nextflow Workflow for Variant Calling
===============================================
*/

nextflow.enable.dsl=2

/*
 * Workflow Input Parameters
 */

params.project = "/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
params.aligned_reads_bam = "${params.project}/Results/Results_Subsample_*/mark_duplicates/duplicates_alignments/*.bam"
params.genome ="${params.project}/Genome/Chr7/chr7.fa"
params.bed = "${params.project}/Metadata/chr7_target.bed"
params.outdir = "${params.project}/Results"


/*
 * Workflow Processes
 */

process LOFREQ_SOMATIC_VC {

    tag {"Variant Calling ${normal_sample} and ${tumor_sample} with Varscan2"}

    container 'quay.io/jeltje/varscan2:latest'

    publishDir "${params.outdir}/${subsample_id}/variant_calling/varscan_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), path(normal_bam), path(normal_bai), val(tumor_sample), path(tumor_bam), path(tumor_bai)
    path fa_file
    path fai_file
    path bed_file

    output:
    path("${subsample_id}*.vcf"), emit: vcf_out, optional: true
    
    script:
    def sample_id_custom = "${normal_sample}_vs_${tumor_sample}"
    """
    varscan2 somatic \\
        -c ${normal_bam} \\
        -t ${tumor_bam} \\
        -q ${sample_id_custom} \\
        -i ${fa_file} \\
        -w ${bed_file} \\
    """
}


/*
 * Workflow
 */

workflow {
    // Definició dels Canals dins del Workflow
    alignment_ch = Channel.fromPath(params.aligned_reads_bam, checkIfExists: true)
        .map { file -> 
            def sub_id = file.parent.parent.parent.name
            def sample_id = file.name.take(file.name.lastIndexOf('.'))
            def type = file.name.take(1) // B, D o E
            return tuple(sub_id, type, sample_id, file, file + ".bai")
        }

    normal_sample_ch = alignment_ch
        .filter { it[1] == 'B' }
        .map { tuple(it[0], it[2], it[3], it[4]) }

    tumor_sample_ch = alignment_ch
        .filter { it[1] != 'B' }
        .map { tuple(it[0], it[2], it[3], it[4]) }

    paired_samples_ch = normal_sample_ch.cross(tumor_sample_ch)
        .map { n, t -> 
            return tuple(n[0], n[1], n[2], n[3], t[1], t[2], t[3])
        }

    // Canals de referència
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true).collect()
    fai_ch    = Channel.fromPath("${params.genome}.fai", checkIfExists: true).collect()
    bed_ch    = Channel.fromPath(params.bed, checkIfExists: true).collect()

    // Execució del Procés
    LOFREQ_SOMATIC_VC(paired_samples_ch, genome_ch, fai_ch, bed_ch)
}