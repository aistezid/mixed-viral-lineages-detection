#!/usr/bin/env python
# coding: utf-8


import csv
import os
from Bio import SeqIO
import gzip

vcf_file = snakemake.input.vcf
reference_file = snakemake.input.ref
output_file = snakemake.output[0]

min_depth = snakemake.params.min_depth
min_alt_freq = snakemake.params.min_alt_freq
sample_name = snakemake.wildcards.sample
reference_seqs = {}


for record in SeqIO.parse(reference_file, "fasta"):
    reference_seqs[record.id] = str(record.seq)

EXPECTED_SEGMENTS = len(reference_seqs) ## number of segments

vcf_data = {}
mixed_positions = []

with gzip.open(vcf_file, "rt") as f:
    for line in f:
        if line.startswith("#"): ## skip header lines
            continue

        fields = line.rstrip().split("\t") ## split vcf line to columns, extract info
        chrom = fields[0]
        pos = int(fields[1])
        ref_nt = fields[3]
        alt_nt = fields[4].split(",")
        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        fmt = dict(zip(format_keys, sample_values)) ## parse the FORMAT column and save it to a dictionary

        ad = list(map(int, fmt["AD"].split(","))) ## split allele depth info to integers
        alleles = [ref_nt] + alt_nt ## get list of all alleles in the site
        allele_counts = dict(zip(alleles, ad)) ## make a dict with alleles and their read coutns
        # print(allele_counts)

        depth = sum(allele_counts.values()) ## calculate total coverage across all alleles
        low_depth = depth < min_depth ## get sites with low depth coverage

        allele_freqs = {}
        high_cov_alleles = {}

        if not low_depth:
            for allele, count in allele_counts.items():
                if count > 0:
                    freq = count / depth ## calculate frequencies of alleles
                    allele_freqs[allele] = freq
                    if freq >= min_alt_freq: ## if frequency is higher or equal than treshold, keep the allele
                        high_cov_alleles[allele] = freq

        if chrom not in vcf_data:
            vcf_data[chrom] = {}

        vcf_data[chrom][pos] = { ## store all calculated info of vcfs in one dict
            "ref": ref_nt,
            "depth": depth,
            "alleles": allele_counts,
            "high": high_cov_alleles,
            "low_depth": low_depth
        }
        # print(vcf_data)

        if len(high_cov_alleles) > 1: ## if there is more than 1 allele of high coverage, save it 
            mixed_positions.append((chrom, pos, ref_nt, depth, high_cov_alleles))

if len(vcf_data) != EXPECTED_SEGMENTS:
    open(output_file, "w").close()
    exit()

## write a report
    output_file = os.path.join(output_path, f"{sample_name}_report.tsv")

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["sample", "segment", "position", "depth", "ref_allele", "allele:count:frequency"]) ## header row

    for segment, pos, ref, depth, alleles in mixed_positions:
        allele_strs = []
        for allele, freq in alleles.items():
            count = int(freq * depth)
            allele_strs.append(f"{allele}:{count}:{freq:.4f}")
            # print(allele_strs)

        writer.writerow([sample_name,segment,pos,depth,ref,",".join(allele_strs)])

