/*
==============================================
    Nextflow Workflow for fastq quality control, alignment, sorting and indexing, marking duplicates
===============================================
*/

nextflow.enable.dsl=2

// Workflow Input Parameters

params.outdir = "results"
params.genome =
params.reads = "Data/Sample_*/*.fastq.gz"


// Workflow Processes
process FASTQC {

    input:

    fastq
    output:

    path "fastqc", emit: fastqc_out
    script:
    """
    fastqc -o fastqc ${input.fastq} -t 2 -q
    """ 

}


process TRIMMING {

}

process BWA_MEM {

}



process SAMTOOLS_SORT {

}


process SAMTOOLS_INDEX {

}

process PICARD_MARKDUPLICATES {

}


// Workflow
workflow {

    main:

    FASTQC(
    
    docker pull biocontainers/fastqc:v0.11.9_cv8

    fastqc_out = FASTQC.out.collect()

    )

    TRIMMING()

    BWA_MEM()

    SAMTOOLS_SORT()

    SAMTOOLS_INDEX()

    PICARD_MARKDUPLICATES()

    publish:
    fastqc_reports = FASTQC.out.collect()
}

output {

    fastqc_reports {
        path "results/fastqc", emit: fastqc_reports
    }
    
}