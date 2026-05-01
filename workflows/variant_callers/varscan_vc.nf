/*
==================================================================================================================================
# VarsCan SOMATIC VARIANT CALLING PIPELINE
==================================================================================================================================
# This script calls somatic SNVs and InDels using the VarsCan variant caller by:
# 1. Calling variants with VarsCan Somatic: variants are called in normal-tumor pairs, defining sample B as normal 
#    and samples D and E as tumoral. Therefore, variant calling is performed for B-D and B-E pairs.
    -n ${normal_bam}: indicate normal sample alignments (B)
    -t ${tumor_bam}: indicate normal sample alignments (D or E)
    -f ${fa_file}: indicate chromosome 7 fasta sequence
    -l ${bed_file}: indicate bed file with Onco Panel sequenced regions on chromosome 7
    --threads 4: number of parallel processes
    -o ${vcf_prefix}: define output vcf file
    --call-indels: call indels with SNVs
# 2. Merging identified SNVs and InDels into a single indexed VCF file

# For each step, a specific Seqera container (compatible with bioconda and arm64/linux) has been used
# ==================================================================================================================================
*/

nextflow.enable.dsl=2

/*
 * Workflow Input Parameters
 */

params.project = "/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
// Path to aligned reads with marked duplicates
params.aligned_reads_bam = "${params.project}/results/Results_Subsample_*/mark_duplicates/duplicates_alignments/*.bam"
// Path to Chromosome 7 reference FASTA and its corresponding index files
params.genome ="${params.project}/genome/Chr7/chr7.fa"
// Path to Onco Panel bed files from Chromosome 7
params.bed = "${params.project}/metadata/bed_files/chr7_target.bed"
// Path to the main Results directory
params.outdir = "${params.project}/results"

/*
 * Workflow Processes
 */

// 1. Calling somatic SNVs and InDels with LoFreq
process samtoolsMpileup {

    tag {"MPileup ${normal_sample} and ${tumor_sample}"}

    container 'community.wave.seqera.io/library/samtools:1.23.1--e8c68bc6da750dc8'
 
    input:
    tuple val(subsample_id), val(normal_sample), path(normal_bam), path(normal_bai), val(tumor_sample), path(tumor_bam), path(tumor_bai)
    path fa_file
    path fai_file
    path bed_file

    output:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path("${tumor_sample}_vs_${normal_sample}.mpileup"), emit: mpileup_ch
    
    script:
    """
    samtools mpileup -f ${fa_file} -l ${bed_file} ${normal_bam} ${tumor_bam} > ${tumor_sample}_vs_${normal_sample}.mpileup
    """
}

// 2. Calling somatic SNVs and InDels with LoFreq
process varscanSomaticVC {

    tag {"Variant Calling ${normal_sample} and ${tumor_sample} with Varscan Somatic"}

    container 'community.wave.seqera.io/library/varscan:2.4.6--f136b57b5c999502'

    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/varscan_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path(mpileup_file)

    output:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path("*.vcf"), emit: vcf_group
    
    script:
    def vcf_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}"
    """
    varscan somatic ${mpileup_file} ${vcf_prefix} \\
        --mpileup 1 \\
        --output-vcf 1 \\
        --min-var-freq 0.01 
    """
}

// 3. Calling somatic SNVs and InDels with LoFreq
process bcftoolsMergeVariants {

    tag {"Filtering ${normal_sample} and ${tumor_sample} variants with bcftools"}

    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'

    publishDir "${params.outdir}/Results_Subsample_${subsample_id.replace('S', '')}/variant_calling/benchmarking_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path(vcf_files)

    output:
    tuple val(subsample_id), path("${subsample_id}_${normal_sample}_vs_${tumor_sample}_varscan_final.vcf.gz*"), emit: varscan_final_vcf
    
    script:
    def merged_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}_varscan_final.vcf.gz"
    """
    for vcf in ${vcf_files}; do
        bcftools view -i 'INFO/SS == 2' \$vcf -Oz -o filtered_\$vcf.gz
        bcftools index -c filtered_\$vcf.gz
    done
    bcftools concat -a filtered_*.gz -Oz -o ${merged_prefix}
    bcftools index ${merged_prefix}
    """
}

/*
 * Workflow
 */

workflow {

    // Channel for aligned .bam .bai files wth marked duplicates for each Subsample 1, 3 and 5
    alignment_ch = Channel.fromPath(params.aligned_reads_bam, checkIfExists: true)
        .map { file -> 
            def sub_id = file.parent.parent.parent.name.replace("Results_Subsample_", "S")
            def duplicates_id = file.name.take(file.name.lastIndexOf('.'))
            def sample_id = duplicates_id.replace("_duplicates", "")
            def type = file.name.take(1)
            return tuple(sub_id, type, sample_id, file, file + ".bai")
        }

    // Channel for Normal samples (B)
    normal_sample_ch = alignment_ch
        .filter { it[1] == 'B' }
        .map { tuple(it[0], it[2], it[3], it[4]) }

    // Channel for Normal samples (D/E)
    tumor_sample_ch = alignment_ch
        .filter { it[1] != 'B' }
        .map { tuple(it[0], it[2], it[3], it[4]) }

    // Pair normal and tumor samples by their Subsample ID
    paired_samples_ch = normal_sample_ch.cross(tumor_sample_ch)
        .map { n, t -> 
            return tuple(n[0], n[1], n[2], n[3], t[1], t[2], t[3])
        }

    // Channel for Chromosome 7 reference FASTA and all index files
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true).collect()
    fai_ch    = Channel.fromPath("${params.genome}.fai", checkIfExists: true).collect()
    // Channel for Onco Panel BED file
    bed_ch    = Channel.fromPath(params.bed, checkIfExists: true).collect()

    // Processes execution
    samtoolsMpileup(paired_samples_ch, genome_ch, fai_ch, bed_ch)
    varscanSomaticVC(samtoolsMpileup.out.mpileup_ch)
    bcftoolsMergeVariants(varscanSomaticVC.out.vcf_group)
}