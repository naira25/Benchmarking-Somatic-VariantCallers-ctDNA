# Benchmarking-Somatic-VariantCallers-ctDNA

This repository contains the scripts and configuration files developed for the **Benchmarking of Somatic Variant Callers in Simulated circulating tumor DNA** Master Thesis. Each script 

The repository is structured as follows:

---

#### Nextflow Configuration Directory (`conf/`)
This directory contains the `nextflow.config` file. This configuration establishes the execution environment (software, containers) and tool-specific parameters required for the Nextflow pipeline's reproducibility.

#### Failed Scripts Directory (`failed_scripts/`)
This directory contains initial approaches and script versions that did not meet the performance or compatibility criteria required for the final pipeline. They are preserved for documentation purposes and to provide context on the problems and adopted optimization process.

#### Processing Scripts Directory (`processing_scripts/`)
This directory contains Bash (`.sh`) scripts for data processing tasks that support or complement the main Nextflow pipelines:
- `genome_hg38/`: Scripts to fetch, decompress and index the *Homo sapiens* (GRCh38) chromosome 7 reference genome.
- `get_fastq_samples/`: Scripts for fetching paired-end raw files from SRA accessions (from samples B, D and E) and managing subsamples files from different sequencing laboratories.
- `truth_set/`: Scripts for standardizing the truth set's variants for comparison and ensuring genomic consistency.
- `variant_calling/`: Scripts to normalize VCF outputs to ensure compatibility between different callers.

#### Nextflow Pipelines Directory (`workflows/`)
This directory contains the core Nextflow pipelines (`.nf`) developed for the main steps of the study:
- `preprocessing/`: Quality control (FastQC), read alignment (BWA-MEM), coordinate sorting, indexing and marking duplicates (GATK/Picard).
- `variant_callers/`: Somatic variant detection using a multi-caller approach, including Mutect2, VarDict, VarScan2 and LoFreq, and an Ensemble approach.
- `benchmarking/`: Automated comparison of caller performance against the SEQC2 Truth Set to calculate sensitivity, precision, and F1-score metrics.

#### .gitignore File
The `.gitignore` file specifies which file extensions, local directories and software are excluded from the GitHub repository. It is configured to keep the repository lightweight by ignoring large genomic data (BAM, FastQ, VCF) while ensuring all source code is tracked.

#### README.md File
This file provides a comprehensive overview of the repository structure, directory contents and necessary documentation to navigate the project.

---

*Note: Further details describing each step, specific tools and command-line parameters are indicated within the comments of each individual script.*
*The scripts and workflows are explicitly cross-referenced throughout the thesis to ensure traceability and easy access to the corresponding source code.*