// FastQC module
process fastqc {
    tag "$sample_id"
    publishDir "$params.outdir/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc.zip"

    script:
    """
    fastqc -o . $reads
    """
}
