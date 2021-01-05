#!/usr/bin/env nextflow


// Parameters
params.outdir = "./"
params.reads = "/data/MSc/2020/MA5112/week_1/raw_data/*_r{1,2}.fastq.gz"

// Initialise Channel
reads_ch = Channel.fromFilePairs(params.reads)


process FastQC {

	label 'FastQC'
	publishDir "${params.outdir}/QC/raw", mode:'copy'

	input:
		tuple val(key), file(reads) from reads_ch

	output:
		file("*.{html,zip}") into fastqc_ch

	script:
	"""
	fastqc -q $reads
	"""
}


process MultiQC {

	label 'MultiQC'
	publishDir "${params.outdir}/QC/raw", mode:'copy'

	input:
		file(htmls) from fastqc_ch.collect()

	output:
		file("*multiqc_report.html") into multiqc_report_ch

	script:
	"""
	multiqc .
	"""
}
