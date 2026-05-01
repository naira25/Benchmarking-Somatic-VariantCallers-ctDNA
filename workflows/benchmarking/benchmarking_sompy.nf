/*
==============================================
    Nextflow Workflow for Somatic Benchmarking
===============================================
*/

nextflow.enable.dsl=2

params.project = "/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
params.genome ="${params.project}/Genome/Chr7/chr7.fa"
params.genome_index = "${params.project}/Genome/Chr7/chr7.fa.fai"
params.input_vcfs = "${params.project}/Results/Results_Subsample_*/variant_calling/benchmarking_vcf/benchmarking_preprocessed/*.vcf.gz"
params.truth = "${params.project}/Metadata/TruthSet/TruthSet_chr7/KnownPositives_chr7_rename.vcf.gz"
params.bed = "${params.project}/Metadata/chr7_target.bed.gz"
params.outdir = "${params.project}/Results"

process vcfStats {
    tag "Stats: ${group_id}"
    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results/${group_id}/vcf_stats", mode: 'copy'

    input:
    tuple val(subsample_id), val(group_id), val(vc_name), path(vcf)

    output:
    tuple val(subsample_id), val(group_id), path("${vc_name}.stats.txt"), emit: stats

    script:
    """
    bcftools stats ${vcf} > ${vc_name}.stats.txt
    """
}

process sompyBenchmark {
    tag "Sompy: ${vc_name}"
    container 'community.wave.seqera.io/library/hap.py:0.3.15--cc9c0286f5a6f629'
    containerOptions "--platform linux/amd64"
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results/${group_id}/sompy_results", mode: 'copy'
    
    input:
    tuple val(subsample_id), val(group_id), val(vc_name), path(test_vcf)
    path truth_vcf
    path truth_tbi
    path fa_file
    path fa_index
    path bed_file

    output:
    tuple val(subsample_id), val(group_id), path("*.stats.csv"), emit: stats
    tuple val(subsample_id), val(group_id), path("*.summary.csv"), emit: summary
    path("*.features.csv"), emit: features
    tuple val(subsample_id), val(group_id), path("${vc_name}_benchmarking.tp.vcf.gz"), emit: tp_vcf, optional: true
    tuple val(subsample_id), val(group_id), path("${vc_name}_benchmarking.fp.vcf.gz"), emit: fp_vcf, optional: true

    script:
    def prefix = "${vc_name}_benchmarking"
    """
    bcftools index -t ${test_vcf}
    export HGREF=${fa_file}
    som.py \\
        ${truth_vcf} \\
        ${test_vcf} \\
        -r ${fa_file} \\
        -R ${bed_file} \\
        -o "${prefix}" \\
        -P \\
        --happy-stats \\
        --feature-table generic \\
        --af-query "AF" \\
        --af-truth "AF"
    """
}

process mergeMetrics {
    tag "Metrics: ${subsample_id}"
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results/${group_id}/merged_parameters", mode: 'copy'
    
    input:
    tuple val(subsample_id), val(group_id), path(summaries)

    output:
    path "parameters_summary.csv"

    script:
    """
    echo "CALLER,SUBSAMPLE,GROUP,TYPE,TRUTH.TOTAL,QUERY.TOTAL,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,METRIC.Frac_NA,METRIC.Ack,QUERY.UNK,FP,TP,FN" > parameters_summary.csv
    for f in ${summaries}; do
        caller=\$(basename \$f | sed 's/_benchmarking.*//')
        # Busquem la línia que conté "SNV" o "Records" (la fila de dades de sompy)
        # i la netegem per si té espais extres
        data=\$(grep -v "Type" \$f | head -n 1)
        echo "\$caller,${subsample_id},${group_id},\$data" >> parameters_summary.csv
    done
    """
}

process multiQCreports {
    tag "MultiQC: ${subsample_id}"
    container 'multiqc/multiqc:latest'
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results/${group_id}/multiQC_reports", mode: 'copy'

    input:
    tuple val(subsample_id), val(group_id), path(summaries), path(stats)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {
    genome_ch = Channel.fromPath(params.genome).collect()
    genome_index_ch = Channel.fromPath(params.genome_index).collect()
    bed_ch = Channel.fromPath(params.bed).collect()
    truth_vcf_ch = Channel.fromPath(params.truth).collect()
    truth_tbi_ch = Channel.fromPath("${params.truth}.tbi").collect()

    test_vcf_ch = Channel.fromPath(params.input_vcfs, checkIfExists: true)
        .map { file -> 
            def sub_id = file.parent.parent.parent.parent.name.replace("Results_Subsample_", "S")
            def group_id = file.name.contains("_D_") ? "B_D" : "B_E"
            def vc_name  = file.baseName.replaceAll(".vcf.gz", "").replaceAll(".vcf", "")
            return tuple(sub_id, group_id, vc_name, file)
        }
    
    vcfStats(test_vcf_ch)
    sompyBenchmark(test_vcf_ch, truth_vcf_ch, truth_tbi_ch, genome_ch, genome_index_ch, bed_ch)

    // MERGE PARAMETERS
    metrics_in = sompyBenchmark.out.summary.groupTuple(by: [0, 1])
    mergeMetrics(metrics_in)
    
    // MULTIQC
    multiqc_in = sompyBenchmark.out.summary
        .groupTuple(by: [0, 1])
        .join(vcfStats.out.stats.groupTuple(by: [0, 1]), by: [0, 1])
    
    multiQCreports(multiqc_in)
}