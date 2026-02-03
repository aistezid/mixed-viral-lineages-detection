# mixed-viral-lineages-detection

This Snakemake pipeline analyses sequencing data to assess whether viral samples contain evidence of multiple viral lineages.

## Input data

The pipeline expects SRA files to be downloaded prior to execution and organised in the samples/ directory. Input files should be placed in separate subdirectories according to sequencing layout (paired-end and single-end reads).

## Processing and filtering criteria

The following filtering thresholds are applied during analysis:

- Reads shorter than 50 nucleotides are discarded

- Reads with a Phred quality score below 20 are excluded

- A minimum sequencing depth of 10 reads is required for allele frequency calculations

- Alternative alleles are considered only if they represent at least 20% of reads at a given position

These thresholds are intended to reduce noise from low-quality reads and low-confidence variants. All thresholds can be adjusted.

## Usage 

Before running the workflow, ensure that:

- Snakemake is installed

- Conda is available

- Input SRA files are organised as described above

To execute the pipeline, run:

snakemake --use-conda