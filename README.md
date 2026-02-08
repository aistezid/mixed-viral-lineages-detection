# mixed-viral-lineages-detection

This Snakemake pipeline analyses sequencing data to assess whether viral samples contain evidence of multiple viral lineages. It is designed for viruses with segmented genomes.

## Requirements

All software dependencies are managed via conda environments defined in the `envs/` directory. Snakemake version 9.14.8 is used.

## Input data

The pipeline expects these inputs: 

- Genome reference FASTA file  located in the `references/` directory

- A comma separated file with SRA accessions including information on whether each sample is paired-end or single-end

## Processing and filtering criteria

Pipeline parameters are defined in `config.yaml` file. The following filtering thresholds are applied by default, but can be modified:

- Reads shorter than 50 nucleotides are discarded

- Reads with a Phred quality score below 20 are excluded

- A minimum sequencing depth of 10 reads is required for allele frequency calculations

- Alternative alleles are considered only if they represent at least 20% of reads at a given position

- Samples must contain all genome segments for sequence assembly

- If more than 5% of nucleotides within a genome segment are ambiguous, the analysis for that sample is stopped

These thresholds are intended to reduce noise from low-quality reads and low-confidence variants. All thresholds can be adjusted.

## Output files

The pipeline produces:

- Tab-separated reports summarizing allele frequencies at highly variable genome positions

- Assembled consensus FASTA sequences in which nucleotides differing from the reference genome are incorporated