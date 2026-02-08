import os
import pandas as pd

REF = "references/wmv6-cms001_028_alco-ref.fasta"
VCF_DIR = "vcf"
samples = pd.read_csv(
	config["sra_list"],
	sep=",",
	header=None,
	names=["SAMPLES", "type"],
	dtype=str,
	keep_default_na=False)

paired_samples = samples[samples["type"] == "paired"]["SAMPLES"].tolist()
single_samples = samples[samples["type"] == "single"]["SAMPLES"].tolist()

SAMPLES=paired_samples + single_samples

rule all:
	input:
		expand("assembled_seq/{sample}_alt.fasta", sample=SAMPLES)

rule download_sra_pe:
	output:
		"samples/paired/{sample}/{sample}.sra"
	shell:
		"prefetch {wildcards.sample}"

rule download_sra_se:
	output:
		"samples/single/{sample}/{sample}.sra"
	shell:
		"prefetch {wildcards.sample}"

rule fastqdump_pe:
	input:
		"samples/paired/{sample}/{sample}.sra"
	output:
		"samples/paired/{sample}/{sample}_1.fastq",
		"samples/paired/{sample}/{sample}_2.fastq"
	conda:
		"envs/fastq.yaml"
	shell:
		"fastq-dump --split-files {input} --outdir samples/paired/{wildcards.sample}"

rule fastqdump_se:
	input:
		"samples/single/{sample}/{sample}.sra"
	output:
		"samples/single/{sample}/{sample}.fastq"
	conda:
		"envs/fastq.yaml"
	shell:
		"fastq-dump {input} --outdir samples/single/{wildcards.sample}"

rule trim_reads_pe:
	input:
		r1="samples/paired/{sample}/{sample}_1.fastq",
		r2="samples/paired/{sample}/{sample}_2.fastq",
	output:
		trimmed_r1="trimmed/paired/{sample}/{sample}_1_trimmed.fastq",
		trimmed_r2="trimmed/paired/{sample}/{sample}_2_trimmed.fastq",
		report_html="trimmed/paired/{sample}/fastq_report.html",
		report_json="trimmed/paired/{sample}/fastq_report.json"
	params:
		min_length=config["min_read_length"],
		quality_threshold=config["quality_threshold"]
	threads: 4
	conda:
		"envs/fastq.yaml"
	shell:
		"""
		mkdir -p trimmed/paired

		fastp \
			--in1 {input.r1} \
			--in2 {input.r2} \
			--out1 {output.trimmed_r1} \
			--out2 {output.trimmed_r2} \
			--thread {threads} \
			--length_required {params.min_length} \
			--qualified_quality_phred {params.quality_threshold} \
			--html {output.report_html} \
			--json {output.report_json} \
			--verbose
		"""

rule trim_reads_se:
	input:
		single="samples/single/{sample}/{sample}.fastq"
	output:
		trimmed_single="trimmed/single/{sample}/{sample}_trimmed.fastq",
		report_html="trimmed/single/{sample}/fastq_report.html",
		report_json="trimmed/single/{sample}/fastq_report.json"
	params:
		min_length=config["min_read_length"],
		quality_threshold=config["quality_threshold"]
	threads: 4
	conda:
		"envs/fastq.yaml"
	shell:
		"""
		mkdir -p trimmed/single

		fastp \
			--in1 {input.single} \
			--out1 {output.trimmed_single} \
			--thread {threads} \
			--length_required {params.min_length} \
			--qualified_quality_phred {params.quality_threshold} \
			--html {output.report_html} \
			--json {output.report_json} \
			--verbose
		"""

rule map_to_reference_pe:
	input:
		r1="trimmed/paired/{sample}/{sample}_1_trimmed.fastq",
		r2="trimmed/paired/{sample}/{sample}_2_trimmed.fastq"
	output:
		bam="mapped/{sample}/{sample}_pe_mapped.bam",
		bai="mapped/{sample}/{sample}_pe_mapped.bam.bai"
	threads: 4
	conda:
		"envs/bwa_samtools.yaml"
	shell:
		"""
		mkdir -p mapped/{wildcards.sample}

		bwa-mem2 index {REF}
		bwa-mem2 mem -t {threads} {REF} {input.r1} {input.r2} | \
		samtools view -bS | \
		samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
		"""

rule map_to_reference_se:
	input:
		single="trimmed/single/{sample}/{sample}_trimmed.fastq"
	output:
		bam="mapped/{sample}/{sample}_se_mapped.bam",
		bai="mapped/{sample}/{sample}_se_mapped.bam.bai"
	threads:4
	conda:
		"envs/bwa_samtools.yaml"
	shell:
		"""
		mkdir -p mapped/{wildcards.sample}

		bwa-mem2 index {REF}
		bwa-mem2 mem -t {threads} {REF} {input.single} | \
		samtools view -bS | \
		samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
		"""

def get_bam(wildcards):
	return [
		f"mapped/{wildcards.sample}/{wildcards.sample}_pe_mapped.bam",
		f"mapped/{wildcards.sample}/{wildcards.sample}_se_mapped.bam",
		]

rule vcf_calling:
	input:
		bam=get_bam
	output:
		vcf="vcf/{sample}/{sample}.vcf.gz"
	conda:
		"envs/bcftools.yaml"
	shell:
		"""
		mkdir -p vcf/{wildcards.sample}

		bcftools mpileup -f {REF} {input.bam} | \
		bcftools call -mv -Oz -o {output.vcf}

		bcftools index {output.vcf}
		"""

rule analyze_vcfs:
	input:
		vcf="vcf/{sample}/{sample}.vcf.gz",
		ref = REF
	output:
		"output/{sample}_report.tsv"
	params:
		min_depth=config["min_depth"],
		min_alt_freq=config["min_alt_freq"]
	conda:
		"envs/vcf_analysis.yaml"
	script:
		"scripts/VCF_frequence_analysis.py"

rule generate_fasta:
	input:
		ref = REF,
		vcf_calc = "output/{sample}_report.tsv"
	output:
		fasta = "assembled_seq/{sample}_alt.fasta"
	conda:
		"envs/vcf_analysis.yaml"
	params:
		min_depth=config["min_depth_assembly"],
		segment_number=config["segment_number"],
		complete_genome=config["complete_genome"],
		max_N_fraction = config["max_N_fraction"]
	script:
		"scripts/snakemake-consensus.py"
