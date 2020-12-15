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

I will walk you through how to construct the processes, leaving you to stitch them together to make a single script.

# Parameters {#params}
In nextflow, parameters are variables that can be passed to the script. You can hard-code them in the script, or assign them a value in the command line.

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

To pass parameters via the command line call:

```bash
nextflow run quality_control.nf --reads "/data/MSc/2020/MA5112/week_1/raw_data/*.fastq.gz"
```

Both produce

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
