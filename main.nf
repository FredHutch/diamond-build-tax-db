#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false

// Containers to use
container__ubuntu = "ubuntu:20.04"
container__prodigal = "quay.io/biocontainers/prodigal:2.6.3--h516909a_2"

// Path to NCBI Taxonomy
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

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

    // Get the NCBI Taxonomy
    get_NCBI_taxonomy(
        file(params.ncbi_taxdump)
    )

    // Extract the CDS from each genome
    prodigal(
        Channel.from(
            file(params.manifest)
        ).splitCsv(
            header: true, 
            sep: ","
        ).map {
            r -> [file(r.uri), r.tax_id]
        }
    )

    // Annotate the gene names with the organism Tax ID
    format_gene_taxid_table(
        prodigal.out
    )

}

// Unpack the NCBI taxonomy and return the nodes.dmp
process get_NCBI_taxonomy {
    container "${container__ubuntu}"
    label "io_limited"
    // errorStrategy "retry"
    
    input:
    file ncbi_taxdump

    output:
    file "nodes.dmp"

"""#!/bin/bash

tar xzvf ${ncbi_taxdump}
"""
}


// Annotation of coding sequences with prodigal
process prodigal {
    container "${container__prodigal}"
    label 'io_limited'
    // errorStrategy "retry"
 
    input:
        tuple file(fasta), val(tax_id)
    
    output:
        tuple file("OUTPUT.faa.gz"), val(tax_id)
    
"""#!/bin/bash

set -e

# Decompress the genome, if it is in gzip format
[[ gzip -t ${fasta} ]] && gunzip -c ${fasta} > INPUT.fasta

# If not, just rename it
[[ gzip -t ${fasta} ]] || cp ${fasta} > INPUT.fasta

prodigal \
    -a OUTPUT.faa \
    -i INPUT.fasta

cat OUTPUT.faa | sed 's/ .*//' | gzip -c > OUTPUT.faa.gz
"""
}

// Make a TSV with the Tax ID for each gene in each genome
process format_gene_taxid_table {
    container "${container__ubuntu}"
    label "io_limited"
    // errorStrategy "retry"
    
    input:
    tuple file(faa_gz), val(tax_id)

    output:
    file "genome_prot2taxid.tsv.gz"

"""#!/bin/bash

gunzip -c "${faa_gz}" | \
    grep '>' | \
    tr -d '>' | \
    while read gene_id; do
        echo \$gene_id \$gene_id ${tax_id} \$gene_id
    done | tr ' ' '\\t' | gzip -c > genome_prot2taxid.tsv.gz
"""
}

// Join together all of those TSVs
process join_gene_taxid_tables {
    container "${container__ubuntu}"
    label "io_limited"
    // errorStrategy "retry"
    
    input:
    file "genome_prot2taxid.*.tsv.gz"

    output:
    file "genome_prot2taxid.tsv.gz"

"""#!/bin/bash

# Write the header
echo "accession	accession.version	taxid	gi" > genome_prot2taxid.tsv

# Add the rows for each genome
for fp in genome_prot2taxid.*.tsv.gz; do
    gunzip -c \$fp >> genome_prot2taxid.tsv
done

gzip genome_prot2taxid.tsv
"""
}
