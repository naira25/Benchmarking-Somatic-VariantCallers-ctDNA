/*
==================================================================================================================================
# BENCHMARKING WITH FILTER SOMATIC VARIANT CALLERS PIPELINE
==================================================================================================================================
*/

nextflow.enable.dsl=2

params.project = "/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
params.genome ="${params.project}/genome/chr7_genome/chr7.fa"
params.genome_index = "${params.project}/genome/chr7_genome/chr7.fa.fai"
params.input_vcfs = "${params.project}/results/Results_Subsample_*/variant_calling/benchmarking_vcf/benchmarking_preprocessed/*.vcf.gz"
params.truth = "${params.project}/metadata/truthset_files/TruthSet_chr7/KnownPositives_chr7_rename.vcf.gz"
params.bed = "${params.project}/metadata/bed_files/chr7_target.bed.gz"
params.outdir = "${params.project}/results"

process filterVCF {
    tag "Filter: ${vc_name}"
    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'
    
    input:
    tuple val(subsample_id), val(group_id), val(vc_name), path(vcf)

    output:
    tuple val(subsample_id), val(group_id), val(vc_name), path("${vc_name}.filtered.vcf.gz"), emit: filtered_vcf

    script:

    if ( vc_name.contains("mutect") ) {
        """
        bcftools filter -i 'FMT/DP >= 1000 && FMT/AF >= 0.005' ${vcf} -Oz -o ${vc_name}.filtered.vcf.gz
        bcftools index -t ${vc_name}.filtered.vcf.gz
        """
    } 
    else if ( vc_name.contains("lofreq") ) {
        """
        bcftools filter -i 'INFO/DP >= 1000 && INFO/AF >= 0.005' ${vcf} -Oz -o ${vc_name}.filtered.vcf.gz
        bcftools index -t ${vc_name}.filtered.vcf.gz
        """
    } 
    else if ( vc_name.contains("varscan") ) {
        """
        bcftools view -h ${vcf} > header.vcf
        bcftools view -H ${vcf} | \\
        awk -F'\\t' 'BEGIN {OFS="\\t"} {
            split(\$10, a, ":"); 
            freq=a[6]; 
            sub(/%/, "", freq); 
            if(a[3] >= 1000 && freq >= 0.5) print \$0
        }' > filtered_data.vcf
        cat header.vcf filtered_data.vcf | bcftools view -Oz -o ${vc_name}.filtered.vcf.gz
        bcftools index -t ${vc_name}.filtered.vcf.gz
        
        rm header.vcf filtered_data.vcf
        """
    }
    else if ( vc_name.contains("vardict") ) {
        """
        bcftools filter -i 'FMT/DP >= 1000 && FMT/AF >= 0.005' ${vcf} -Oz -o ${vc_name}.filtered.vcf.gz
        bcftools index -t ${vc_name}.filtered.vcf.gz
        """
    }
    else if ( vc_name.contains("ensemble") ) {
        """
        bcftools filter -i 'FMT/DP >= 1000 && INFO/AF >= 0.005' ${vcf} -Oz -o ${vc_name}.filtered.vcf.gz
        bcftools index -t ${vc_name}.filtered.vcf.gz
        """
    }
}

process vcfStats {
    tag "Stats: ${group_id}"
    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results_filters/${group_id}/vcf_stats", mode: 'copy'

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
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results_filters/${group_id}/sompy_results", mode: 'copy'
    
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
    def af_field = (vc_name.contains("varscan")) ? "FREQ" : "AF"
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
        --af-query "${af_field}" \\
        --af-truth "AF"
    """
}

process mergeMetrics {
    tag "Metrics: ${subsample_id}"
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results_filters/${group_id}/merged_parameters", mode: 'copy'
    
    input:
    tuple val(subsample_id), val(group_id), path(summary_files)

    output:
    path "parameters_summary.csv"

    script:
    """
    echo "CALLER,SUBSAMPLE,GROUP,TRUTH_TOTAL,QUERY_TOTAL,TP,FP,FN,RECALL,PRECISION,F1" > parameters_summary.csv
    for f in ${summary_files}; do
        caller=\$(basename \$f | cut -d'_' -f9)
        line=\$(grep ",PASS," \$f | awk -F',' '{print \$4","\$7","\$5","\$8","\$6","\$11","\$12","\$14}')
        echo "\$caller,${subsample_id},${group_id},\$line" >> parameters_summary.csv
    done
    """
}

process multiQCreports {
    tag "MultiQC: ${subsample_id}"
    container 'multiqc/multiqc:latest'
    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_results_filters/${group_id}/multiQC_reports", mode: 'copy'

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
            def group_id = (file.name.contains("_D_") || file.name.contains("_vs_D_")) ? "B_D" : "B_E"
            def vc_name  = file.baseName.replaceAll(".vcf.gz", "").replaceAll(".vcf", "")
            return tuple(sub_id, group_id, vc_name, file)
        }
    
    vcfStats(test_vcf_ch)
    filterVCF(test_vcf_ch)
    sompyBenchmark(filterVCF.out.filtered_vcf, truth_vcf_ch, truth_tbi_ch, genome_ch, genome_index_ch, bed_ch)
    metrics_in = sompyBenchmark.out.summary.groupTuple(by: [0, 1])
    mergeMetrics(metrics_in)
    multiqc_in = sompyBenchmark.out.summary
        .groupTuple(by: [0, 1])
        .join(vcfStats.out.stats.groupTuple(by: [0, 1]), by: [0, 1])
    multiQCreports(multiqc_in)
}