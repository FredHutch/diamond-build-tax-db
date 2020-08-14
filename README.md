# diamond-build-tax-db
Build DIAMOND database with taxonomic annotations
Author: Sam Minot, Ph.D.

Nextflow utility to build a DIAMOND database using a set of genomes
which have been annotated with an NCBI Tax ID.

### Manifest

All of the input genomes (in FASTA format) must be listed in a manifest,
which will be used as the input to this utility. That CSV must contain
at least two named columns:

- uri: The location of the genomes, either on the local filesystem or in object storage.
- tax_id: The NCBI Tax ID. Make sure that this matches the current taxonomy.

As with any file which may be made with Excel, some of the primary failure
cases for users will come from the way that this file is formatted. Make sure that
the column name matches exactly _after_ exporting to CSV from Excel. It's also
always a good idea to run `dos2unix` to make sure that all carriage returns are
corrected.

### Building the Database

To build the database using a manifest file `manifest.csv` and writing the output
to `my_output_folder/my_database.dmnd`, simply run the command:

```#!/bin/bash

NXF_VER=20.07.0 \
nextflow \
    run \
    FredHutch/diamond-build-tax-db \
    --manifest manifest.csv \
    --output_folder my_output_folder \
    --output_prefix my_database \
    -with-report \
    -with-trace \
    -resume \
    -latest

```

Depending on your compute infrastructure, you may require additional parameters,
most notably indicating your personal `nextflow.config` file with `-c ~/nextflow.config`, e.g.

### Get in Touch

If you have any questions, comments, or concerns, please get in touch with the
maintainer of this repository: Samuel Minot, Ph.D. (sminot at fredhutch dot org).
