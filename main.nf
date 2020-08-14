#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/diamond-build-tax-db <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing genome FASTAs (uri) and NCBI Tax IDs (tax_id)
      --output_folder       Folder to write database to
      --output_prefix       Prefix to use for database

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || !params.manifest || !params.output_folder || !params.output_prefix){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

}