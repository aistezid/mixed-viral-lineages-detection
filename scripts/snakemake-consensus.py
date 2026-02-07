#!/usr/bin/env python
# coding: utf-8


from Bio import SeqIO
import os

min_depth_assembly = snakemake.params.min_depth
expected_segments = snakemake.params.segment_number
require_complete = snakemake.params.complete_genome
max_N_fraction = snakemake.params.max_N_fraction

reference_file = snakemake.input.ref
vcf_output = snakemake.input.vcf_calc
fasta_output = snakemake.output.fasta

reference_seqs = {}
ref_index = {}
for record in SeqIO.parse(reference_file, "fasta"): ## make a dict which stores reference sequence nucleotide positions (per segment)
    segment = record.id.split("|")[1]
    ref_index[segment] = {}
    position = 0
    for nt in record.seq:
        position += 1
        ref_index[segment][position] = nt
    # print(ref_index)

vcf_data = {}
with open(vcf_output, 'r') as f: ## make a dict which stores alternative allele info: position, what is the alternative dict
    for line in f:
        if line.startswith("sample"):
            continue
        segment = line.split("\t")[1].split("|")[1]
        position, depth, ref_allele, alt_allele_info = int(line.split("\t")[2]), line.split("\t")[3], line.split("\t")[4], line.split("\t")[5]
        alt_allele = alt_allele_info.split(",")[1].split(":")[0]
        if segment not in vcf_data:
            vcf_data[segment] = {}
        vcf_data[segment][position] = {"ref_allele": ref_allele, "alt_allele": alt_allele}
    # print(vcf_data)

total_segments = len(ref_index)  # number of segments in the sample
# Decide if this sample is "complete enough" to continue
if require_complete and total_segments < expected_segments:
    continue_assembly = False
else:
    continue_assembly = True

new_seq = {} ## make a dict for storing new sequences
sra_num = vcf_output.split("/")[-1].split("_")[0]

if continue_assembly:
    for segment in ref_index:
        seq = []
        for position in sorted(ref_index[segment]):
            nt = ref_index[segment][position]
            if segment not in vcf_data or position not in vcf_data[segment]: ## if nucleotide position is not in alternative allele dict, then take reference sequence nucleotide
                seq.append(nt)
            else: ## if nt position is in alt allele dict, take alternative nucleotide
                variant = vcf_data[segment][position]
                if variant.get("depth", 0) < 10:
                    seq.append("N")  ## if coverage is low, insert N
                else:
                    seq.append(variant["alt_allele"]) 
        
        N_count_seg = seq.count("N")
        seg_length = len(seq)
        N_fraction = N_count_seg / seg_length

        if N_fraction > max_N_fraction: ## if fraction of ambiguous nucleotides exceeds set fraction, skip the segment
            print(f"Sample {sra_num}, segment {segment} skipped: {N_fraction} N fraction in the segment exceeds max_N_fraction = {max_N_fraction}")
            continue  
        
        header = f"{sra_num}", "|", f"{segment}"
        header_str = "".join(header)
        new_seq[segment] = "".join(seq)
else:
    print(f"Sequence assembly skipped for sample {sra_num} due to incomplete genome")

with open(fasta_output, "w") as f: ## save new sequences with alternative variant positions
    for segment, seq in new_seq.items():
        header = f"{sra_num}|{segment}" ## header format: SRA|segment
        f.write(f">{header}\n{seq}\n")
