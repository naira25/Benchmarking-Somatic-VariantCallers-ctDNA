/*
==================================================================================================================================
# Mutect2 SOMATIC VARIANT CALLING PIPELINE
==================================================================================================================================
# This script calls somatic SNVs and InDels using the GATK4 Mutect2 variant caller by:
# 1. Calling variants with Mutect2: variants are called in normal-tumor pairs, defining sample B as normal 
#    and samples D and E as tumoral. Therefore, variant calling is performed for B-D and B-E pairs.
    -R ${fa_file}: indicates path to chromosome 7 reference FASTA sequence
    -I ${tumor_bam}: indicates path to tumor sample alignments (D or E)
    -I ${normal_bam}: indicates path to normal sample alignments (B)
    -normal: defines the normal sample name for somatic modeling
    -L ${bed_file}: indicates path to BED file defining Onco Panel target regions on chromosome 7
# 2. Filtering and Selection: variants called by Mutect2 are processed through:
#    - FilterMutectCalls: applies technical filters to the raw calls.
#    - SelectVariants: extracts high-confidence SNVs and InDels that passed all filters (--exclude-filtered).
# 3. The final VCF file is renamed and indexed for downstream benchmarking.

# For each step, a specific Seqera container (compatible with bioconda and arm64/linux) has been used
# ==================================================================================================================================
*/

nextflow.enable.dsl=2

/*
 * Workflow Input Parameters
 */

// Path to the main project directory
params.project = "/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"
// Path to aligned reads with marked duplicates
params.aligned_reads_bam = "${params.project}/results/Results_Subsample_*/mark_duplicates/duplicates_alignments/*.bam"
// Path to Chromosome 7 reference FASTA and its corresponding index files
params.genome ="${params.project}/genome/chr7_genome/chr7.fa"
// Path to Onco Panel bed files from Chromosome 7
params.bed = "${params.project}/metadata/bed_files/chr7_target.bed"
// Path to the main Results directory
params.outdir = "${params.project}/results"


/*
 * Workflow Processes
 */
// 1. Calling somatic SNVs and InDels with Mutect2
process mutectSomaticVC {

    tag {"Variant Calling ${normal_sample} and ${tumor_sample} with Mutect2 Somatic"}
    container 'community.wave.seqera.io/library/gatk4:4.6.2.0--eb4eddc44dc7fb63'
    publishDir "${params.outdir}/Results_Subsample_*/variant_calling/mutect_vcf", mode: 'copy'
 
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
    // Define the output directory name
    def vcf_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}"
    """
    # Extract the normal sample name from the BAM header using GATK GetSampleName
    gatk GetSampleName -I ${normal_bam} -O normal_name.txt
    REAL_NORMAL=\$(cat normal_name.txt)
    
    # Run Mutect2 for somatic variant calling
    gatk --java-options "-Xmx3g" Mutect2 \\
        -R ${fa_file} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal "\$REAL_NORMAL" \\
        -L ${bed_file} \\
        -O ${vcf_prefix}.vcf.gz \\
        --native-pair-hmm-threads 2

    # Filter and Select Variants
    gatk --java-options "-Xmx3g"  FilterMutectCalls \\
        -R ${fa_file} \\
        -V ${vcf_prefix}.vcf.gz \\
        --stats ${vcf_prefix}.vcf.gz.stats \\
        -O ${vcf_prefix}_filtered.vcf.gz

    # Select high-confidence SNVs and InDels that passed all filters
    gatk --java-options "-Xmx3g" SelectVariants \\
        -R ${fa_file} \\
        -V ${vcf_prefix}_filtered.vcf.gz \\
        --select-type-to-include SNP \\
        --select-type-to-include INDEL \\
        --exclude-filtered true \\
        -O ${vcf_prefix}_final.vcf.gz
    """
}

// 2. Indexing the final VCF files with bcftools
process bcftoolsIndexVariants {

    tag {"Indexing ${normal_sample} and ${tumor_sample} variants with bcftools"}
    container 'community.wave.seqera.io/library/bcftools:1.23.1--15a480527db1d585'
    publishDir "${params.outdir}/Results_Subsample_*/variant_calling/benchmarking_vcf", mode: 'copy'
 
    input:
    tuple val(subsample_id), val(normal_sample), val(tumor_sample), path(vcf_file)

    output:
    tuple val(subsample_id), path("${subsample_id}_${normal_sample}_vs_${tumor_sample}_mutect_final.vcf.gz*"), emit: mutect_final_indexed
    
    script:
    // Define the output file name for the indexed VCF
    def indexed_prefix = "${subsample_id}_${normal_sample}_vs_${tumor_sample}_mutect_final.vcf.gz"
    """
    # Rename the final VCF file to include the caller name and index it with bcftools
    mv ${vcf_file} ${indexed_prefix}
    bcftools index -c ${indexed_prefix}
    """
}

/*
 * Workflow
 */

workflow {
    // Channel for aligned .bam/.bai files wth marked duplicates for each Subsample 1, 3 and 5
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
    dict_ch   = Channel.fromPath("${params.genome.replace('.fa', '.dict')}", checkIfExists: true).collect()
    // Channel for Onco Panel BED file
    bed_ch    = Channel.fromPath(params.bed, checkIfExists: true).collect()

    // Processes execution
    mutectSomaticVC(paired_samples_ch, genome_ch, fai_ch, dict_ch, bed_ch)
    bcftoolsIndexVariants(mutectSomaticVC.out.vcf_group)
}