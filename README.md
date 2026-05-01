# Benchmarking-Somatic-VariantCallers-synthetic-ctDNA

This repository contains the scripts and configuration files used in the Master's Thesis: **"Benchmarking of Somatic Variant Callers in Synthetic Circulating DNA"**.

The repository is structured as follows:

---

#### Nextflow Configuration Directory (*conf*)
This directory contains the `nextflow.config` file. This configuration establishes the execution environment (profiles, containers) and tool-specific parameters required for the pipeline's reproducibility across different computing platforms.

#### Failed Scripts Directory (*failed_scripts*)
This directory contains initial approaches and script versions that did not meet the performance or compatibility criteria required for the final pipeline. They are preserved for documentation purposes and to provide context on the troubleshooting and optimization process.

#### Processing Scripts Directory (*processing_scripts*)
This directory contains Bash (`.sh`) scripts for data processing tasks that support or complement the main Nextflow pipelines:
- **Genome Data Processing:** Scripts to fetch, decompress, and index the Chromosome 7 reference genome (GRCh38).
- **Raw FastQ Processing:** Fetching paired-end raw FastQ files from SRA accessions (samples B, D, and E) and managing subsampled files from different sequencing laboratories.
- **Truth Set Processing:** Standardizing the "Gold Standard" variants for comparison and ensuring genomic coordinate consistency.
- **Variant Calling Processing:** Scripts to normalize VCF outputs (standardizing headers and sample names) to ensure compatibility between different callers.
- **Benchmarking Processing:** Scripts to implement the ensemble approach and generate a consensus summary of called variants.

#### Workflows Scripts Directory (*workflows*)
This directory contains the core Nextflow Workflows (`.nf`) developed for the main stages of the study:
- **Data Preprocessing:** Quality control (FastQC), read alignment (BWA-MEM), coordinate sorting, indexing, and marking duplicates (GATK/Picard).
- **Variant Calling:** Somatic variant detection using a multi-caller approach, including MuTect2, VarDict, VarScan, and LoFreq.
- **Benchmarking:** Automated comparison of caller performance against the SEQC2 Truth Set to calculate sensitivity, precision, and F1-score metrics.

#### .gitignore File
The `.gitignore` file specifies which file extensions, local directories, and software are excluded from the repository. It is configured to keep the repository lightweight by ignoring large genomic data (BAM, FastQ, VCF) while ensuring all source code is tracked.

#### README.md File
This file provides a comprehensive overview of the repository structure, directory contents, and the necessary documentation to navigate the project.

---

*Note: Further details describing each step, specific tools and command-line parameters are indicated within the comments of each individual script.*