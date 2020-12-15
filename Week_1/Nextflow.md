---
title:
layout: page
permalink: /Week_1/Nextflow
---

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/week1/nextflow.png" width="100%" height="100%"/>
</center>

We will use `nextflow` to construct a script capable of performing read QC, adapter/read filtering and subsequent QC on the trimmed reads.

Jump to:
- [Parameters](#params)
- [Channels](#channels)
- [Exercise](#exercise)

I will show you how to construct the processes, leaving you to stitch them together to make a single script.

# Parameters {#params}
In nextflow, parameters are variables that can be passed to the script. You can hard-code them in the script, or pass them to the script via the command line.

```nextflow
#!/usr/bin/env nextflow

// Parameters
params.reads = "/data/MSc/2020/MA5112/week_1/raw_data/*_r{1,2}.fastq.gz"

// Initialise Channel
reads_ch = Channel.fromFilePairs(params.reads)
reads_ch.view()
```

Run the hard-coded paramters below by calling:

```bash
nextflow run quality_control.nf
```

Pass parameters via the command line by calling:

```bash
nextflow run quality_control.nf --reads "/data/MSc/2020/MA5112/week_1/raw_data/*_r{1,2}.fastq.gz"
```

Both produce the same output

```
N E X T F L O W  ~  version 20.10.0
Launching `quality_control.nf` [berserk_gutenberg] - revision: eace5429d7
[hcc1395_normal_rep1, [/data/MSc/2020/MA5112/week_1/raw_data/hcc1395_normal_rep1_r1.fastq.gz, /data/MSc/2020/MA5112/week_1/raw_data/hcc1395_normal_rep1_r2.fastq.gz]]
[hcc1395_normal_rep2, [/data/MSc/2020/MA5112/week_1/raw_data/hcc1395_normal_rep2_r1.fastq.gz, /data/MSc/2020/MA5112/week_1/raw_data/hcc1395_normal_rep2_r2.fastq.gz]]
[hcc1395_normal_rep3, [/data/MSc/2020/MA5112/week_1/raw_data/hcc1395_normal_rep3_r1.fastq.gz, /data/MSc/2020/MA5112/week_1/raw_data/hcc1395_normal_rep3_r2.fastq.gz]]
[hcc1395_tumor_rep1, [/data/MSc/2020/MA5112/week_1/raw_data/hcc1395_tumor_rep1_r1.fastq.gz, /data/MSc/2020/MA5112/week_1/raw_data/hcc1395_tumor_rep1_r2.fastq.gz]]
[hcc1395_tumor_rep2, [/data/MSc/2020/MA5112/week_1/raw_data/hcc1395_tumor_rep2_r1.fastq.gz, /data/MSc/2020/MA5112/week_1/raw_data/hcc1395_tumor_rep2_r2.fastq.gz]]
[hcc1395_tumor_rep3, [/data/MSc/2020/MA5112/week_1/raw_data/hcc1395_tumor_rep3_r1.fastq.gz, /data/MSc/2020/MA5112/week_1/raw_data/hcc1395_tumor_rep3_r2.fastq.gz]]

```

# Channels {#channel}
In nextflow files are passed to each process via channels.

In the previous example we created the channel `reads_ch` by specifying `Channel.fromFilePairs()`. The channel now contains each file matching the glob pattern `*.fastq.gz` in the directory given by `params.reads`.

We will use the `read_ch` to provide the `FastQC` process the raw data, collect the files prodcued from FastQC and pass them to MultiQC.

*N.B: Channels can only be used once!*

```nextflow
#!/usr/bin/env nextflow


// Parameters
params.outdir = "./"
params.reads = "/data/MSc/2020/MA5112/week_1/raw_data/*_r{1,2}.fastq.gz"

// Initialise Channel
reads_ch = Channel.fromFilePairs(params.reads)

process FastQC {

	label 'FastQC'
	publishDir "$params.outdir/QC/raw", mode:'copy'

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

```

To run the script on Lugh we will add additional parameters to the command line:

```bash
nextflow -bg -q run quality_control.nf -with-singularity container/week1.img
```

1. `-bg`: Run the script in the background. Nextflow will parse the Configuration file located at `~/.nextflow/config` to determine the resource usage for SLURM.
2. `-q`: Quiet, do not print stdout to console.
3. `-with-singularity`: Use container to execute the script. Provide the path to the container.

Run this line of code yourself, and type `squeue -u mscstudent` on Lugh to view the submitted jobs.

When the job is finished, the reports will be under `PATH_TO_YOUR_DIRECTORY/QC/raw`.

# Exercise {#exercise}
***
I have shown you how to intitialise a raw read channel for downstream use with FastQC and MultiQC.

In the exercise below I have started a nextflow script to read in the raw reads and perform adapter trimming and read filtering using bbduk.sh. Finish the script by using the trimmed reads channel, passing it to a FastQC and MultiQC processes to inspect the statistics of the filtered reads.

```nextflow
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


	output:


	script:
	"""
	fastqc -q $reads
	"""
}

process MultiQC {

	label 'MultiQC'
	publishDir "${params.outdir}/QC/trimmed", mode:'copy'

	input:



	output:


	script:
	"""
	multiqc .
	"""
}
```

### To Do:
1. Fill in the inputs and outputs for the process FastQC.
2. Fill in the inputs and outputs for the process MultiQC. Use the html reports **and** statistics text file from bbduk as inputs for this process. (*Hint use .collect() on both inputs*).
3. Save the script and run using `nextflow -bg -q trim_qc.nf -with-singularity container/week1.img`.

***

Congratulations on making it through the week 1 tutorial. Please do not hesistate to contact me if you have any questions! [bdigby1@nuigalway.ie](b.digby1@nuigalway.ie)

Files used for the tutorial are available at the following [link](https://github.com/BarryDigby/barrydigby.github.io/tree/master/Week_1)

> Barry
