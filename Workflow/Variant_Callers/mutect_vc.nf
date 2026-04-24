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

process mutectSomaticVC {

    tag {"Variant Calling ${normal_sample} and ${tumor_sample} with Mutectc Somatic"}

    container 'community.wave.seqera.io/library/gatk4:4.6.2.0--eb4eddc44dc7fb63'

    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/mutect_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), path(normal_bam), path(normal_bai), val(tumor_sample), path(tumor_bam), path(tumor_bai)
    path fa_file
    path fai_file
    path dict_file
    path bed_file

    output:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path("*_final.vcf.gz"), emit: vcf_group
    path("*.stats"), emit: mutect_stats

    script:
    def vcf_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}"
    """
    gatk GetSampleName -I ${normal_bam} -O normal_name.txt
    REAL_NORMAL=\$(cat normal_name.txt)

    gatk --java-options "-Xmx3g" Mutect2 \\
        -R ${fa_file} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal "\$REAL_NORMAL" \\
        -L ${bed_file} \\
        -O ${vcf_prefix}.vcf.gz \\
        --native-pair-hmm-threads 2

    gatk --java-options "-Xmx3g"  FilterMutectCalls \\
        -R ${fa_file} \\
        -V ${vcf_prefix}.vcf.gz \\
        --stats ${vcf_prefix}.vcf.gz.stats \\
        -O ${vcf_prefix}_filtered.vcf.gz

    gatk --java-options "-Xmx3g" SelectVariants \\
        -R ${fa_file} \\
        -V ${vcf_prefix}_filtered.vcf.gz \\
        --select-type-to-include SNP \\
        --select-type-to-include INDEL \\
        --exclude-filtered true \\
        -O ${vcf_prefix}_final.vcf.gz
    """
}

process bcftoolsIndexVariants {

    tag {"Indexing ${normal_sample} and ${tumor_sample} variants with bcftools"}

    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'

    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path(vcf_file)

    output:
    tuple val(subsample_id), path("${subsample_id}_${normal_sample}_vs_${tumor_sample}_mutect_final.vcf.gz*"), emit: mutect_final_indexed
    
    script:
    def indexed_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}_mutect_final.vcf.gz"
    """
    mv ${vcf_file} ${indexed_prefix}
    bcftools index -c ${indexed_prefix}
    """
}

/*
 * Workflow
 */

workflow {
    alignment_ch = Channel.fromPath(params.aligned_reads_bam, checkIfExists: true)
        .map { file -> 
            def sub_id = file.parent.parent.parent.name.replace("Results_Subsample_", "S")
            def duplicates_id = file.name.take(file.name.lastIndexOf('.'))
            def sample_id = duplicates_id.replace("_duplicates", "")
            def type = file.name.take(1)
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
    dict_ch   = Channel.fromPath("${params.genome.replace('.fa', '.dict')}", checkIfExists: true).collect()
    bed_ch    = Channel.fromPath(params.bed, checkIfExists: true).collect()

    mutectSomaticVC(paired_samples_ch, genome_ch, fai_ch, dict_ch, bed_ch)
    bcftoolsIndexVariants(mutectSomaticVC.out.vcf_group)
}