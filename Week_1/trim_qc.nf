#!/usr/bin/env nextflow


// Parameters
params.outdir = "./"
params.reads = "/data/MSc/2020/MA5112/week_1/raw_data/*_r{1,2}.fastq.gz"
params.adapters = "/data/MSc/2020/MA5112/week_1/assets/adapters.fa"

// Initialise Channel
reads_ch = Channel.fromFilePairs(params.reads)

process Trim {

	label 'BBDUK'
	publishDir "${params.outdir}/trimmed_reads", mode:'copy', pattern: "*.fq.gz"
	publishDir "${params.outdir}/QC/trimmed", mode:'copy', pattern: "*.stats.txt"

	input:
		tuple val(key), file(reads) from reads_ch
		path(adapters) from params.adapters

	output:
		tuple val(key), file("*.fq.gz") into trimmed_reads_ch
		file("*.stats.txt") into trimmed_stats_ch

	script:
	"""
	bbduk.sh -Xmx4g \
	in1=${reads[0]} \
	in2=${reads[1]} \
	out1=${key}_r1.fq.gz \
	out2=${key}_r2.fq.gz \
	ref=$adapters \
	minlen=30 \
	ktrim=r \
	k=12 \
	qtrim=r \
	trimq=20 \
	stats=${key}.stats.txt
	"""
}

process FastQC {

	label 'FastQC'
	publishDir "${params.outdir}/QC/trimmed", mode:'copy'

	input:
		tuple val(key), file(reads) from trimmed_reads_ch

	output:
		file("*.{html,zip}") into fastqc_ch

	script:
	"""
	fastqc -q $reads
	"""
}

process MultiQC {

	label 'MultiQC'
	publishDir "${params.outdir}/QC/trimmed", mode:'copy'

	input:
		file(htmls) from fastqc_ch.collect()
		file(stats) from trimmed_stats_ch.collect()
	
	output:
		file("*multiqc_report.html") into multiqc_report_ch

	script:
	"""
	multiqc .
	"""
}
