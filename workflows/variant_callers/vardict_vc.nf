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

process vardictSomaticVC {

    tag {"Variant Calling ${normal_sample} and ${tumor_sample} with VarDict Somatic"}

    container 'community.wave.seqera.io/library/vardict-java:1.8.3--675bbcc5303f55d1'

    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/vardict_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), path(normal_bam), path(normal_bai), val(tumor_sample), path(tumor_bam), path(tumor_bai)
    path fa_file
    path fai_file
    path bed_file

    output:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path("*_raw.vcf"), emit: vcf_group
    path("*.log"), emit: log_out
    
    script:
    def vcf_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}"
    def af_threshold = 0.01
    
    """
    mkdir -p getopt_temp
    Rscript -e "install.packages('getopt', repos='https://cloud.r-project.org', lib='getopt_temp')"
    export R_LIBS_USER=\$(pwd)/getopt_temp
    
    vardict-java \\
        -G ${fa_file} \\
        -f ${af_threshold} \\
        -N "${tumor_sample}" \\
        -b "${tumor_bam}|${normal_bam}" \\
        -c 1 -S 2 -E 3 ${bed_file} \\
    | testsomatic.R \\
    | var2vcf_paired.pl -N "${tumor_sample}|${normal_sample}" -f ${af_threshold} \\
    > ${vcf_prefix}_raw.vcf 2> ${vcf_prefix}.log

    rm -rf getopt_temp
    """
}

process bcftoolsIndexVariants {

    tag {"Indexing ${normal_sample} and ${tumor_sample} variants with bcftools"}

    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'

    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path(vcf_file)

    output:
    tuple val(subsample_id), path("*.vcf.gz*"), emit: vardict_final_indexed
    
    script:
    def indexed_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}_vardict_final.vcf.gz"
    """
    bcftools filter -i 'STATUS ~ "Somatic" && FILTER == "PASS"' ${vcf_file} -Oz -o ${indexed_prefix}
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

    genome_ch = Channel.fromPath(params.genome, checkIfExists: true).collect()
    fai_ch    = Channel.fromPath("${params.genome}.fai", checkIfExists: true).collect()
    bed_ch    = Channel.fromPath(params.bed, checkIfExists: true).collect()

    // Execució del Procés
    vardictSomaticVC(paired_samples_ch, genome_ch, fai_ch, bed_ch)
    bcftoolsIndexVariants(vardictSomaticVC.out.vcf_group)
}