!/usr/bin/env nextflow

nextflow enable.dsl=2

*/
=======================================================================================
    Pipeline for the Benchmarking of Somatic Variant Callers Master Thesis
=======================================================================================
*/

*/
=======================================================================================
    Define the input fastq files corresponding to Onco Panel Sequencing of ctDNA samples B, D and E
=======================================================================================
*/

params.input = "/workspaces/Benchmarking-Somatic-VariantCallers-ctDNA/Samples/Sample_*/fastq_data/fastq/*.fastq.gz"

*/
=======================================================================================
    Define input fastq files through channels 
=======================================================================================
*/

input_channel = Channel.fromPath(params.input)

*/
=======================================================================================
    Main Workflow 
=======================================================================================
*/

workflow {

// Define the main workflow process
    main_process(input_channel)
}