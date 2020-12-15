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

I will show you how to construct the processes, leaving you to stitch them together to make a single script.

# Parameters {#params}
In nextflow, parameters are variables that can be passed to the script. You can hard-code them in the script, or pass them to the script via the command line.

```nextflow
#!/usr/bin/env nextflow

// Parameters
params.reads = "/data/MSc/2020/MA5112/week_1/raw_data/*.fastq.gz"

// Initialise Channel
reads_ch = Channel.fromPath(params.reads)
reads_ch.view()
```

Run the hard-coded paramters below by calling:

```bash
nextflow run quality_control.nf
```

Pass parameters via the command line by calling:

```bash
nextflow run quality_control.nf --reads "/data/MSc/2020/MA5112/week_1/raw_data/*.fastq.gz"
```

Both produce the same output

```
N E X T F L O W  ~  version 20.10.0
Launching `quality_control.nf` [shrivelled_pare] - revision: 8eedf94aa7
/data/MSc/2020/MA5112/week_1/raw_data/A375_1.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/A375_2.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/A375_3.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/A549_1.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/A549_2.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/A549_3.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/ctrl_1.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/ctrl_2.fastq.gz
/data/MSc/2020/MA5112/week_1/raw_data/ctrl_3.fastq.gz
```

# Channels {#channel}
In nextflow files are passed to each process via channels.

In the previous example we created the channel `reads_ch` by specifying `Channel.fromPath()`. The channel now contains each file matching the glob pattern `*.fastq.gz` in the directory given by `params.reads`.

We will use the `read_ch` to provide the `FastQC` process the raw data.

*N.B: Channels can only be used once!*

```nextflow
#!/usr/bin/env nextflow


// Parameters
params.outdir = "./"
params.reads = "/data/MSc/2020/MA5112/week_1/raw_data/*.fastq.gz"

// Initialise Channel
reads_ch = Channel.fromPath(params.reads)

process FastQC {

	label 'FastQC'
	publishDir "$params.outdir/QC/raw", mode:'copy'

	input:
		file(reads) from reads_ch

	output:
		file("*.{html,zip}") into fastqc_ch

	script:
	"""
	fastqc -q $reads
	"""
}
```

To run the script on Lugh we will add additional parameters to the command line:

```bash
nextflow -bg -q run quality_control.nf -with-singularity containers/week1.img
```

1. `-bg`: Run the script in the background. Nextflow will parse the Configuration file located at `~/.nextflow/config` to determine the resource usage for SLURM.
2. `-q`: Quiet, do not print stdout to console.
3. `-with-singularity`: Use container to execute the script. Provide the path to the container.

Run this line of code yourself, and type `squeue -u mscstudent` on Lugh to view the submitted jobs.
