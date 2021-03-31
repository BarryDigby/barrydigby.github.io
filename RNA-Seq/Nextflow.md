---
title: Nextflow
layout: page
permalink: /RNA-Seq/Nextflow
---

# Introduction
There will be no exercise for you to complete in nextflow this week. Instead, I have designed a simple nextflow script that can perform both single-end and paired-end analysis with kallisto by using `if/else` statements.

Furthermore, the script is capable of downloading the reference cDNA index file and performing `kallisto index` if neither file is not provided.

## Single end data
Before we walk thorugh the script, you need to know how to handle single-end fastq files with nextflow. In the example below, we specify the paths to the fastq files using `fromPath()`. This will return a channel with the paths to the fastq files.

Save the below script and run it:

```bash
#!/usr/bin/env nextflow

params.input = "/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/*.fastq.gz"

Channel
      .fromPath(params.input)
      .set{ch_reads}

ch_reads.view()
```

```
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/ctrl_2.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/ctrl_1.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/ctrl_3.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A375_1.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A375_2.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A375_3.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A549_1.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A549_2.fastq.gz
/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A549_3.fastq.gz
```

The channel outputs file paths. We want it to mimick `fromFilePairs()` and create a tuple with a key for each file path, so we can name the output directories when running `kallisto quant` for single end data.

Save the below script and run it:

```bash
#!/usr/bin/env nextflow

params.input = "/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/*.fastq.gz"

Channel
      .fromPath(params.input)
      .map{ file -> [file.simpleName, file]}
      .set{ch_reads}

ch_reads.view()
```

```
[ctrl_2, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/ctrl_2.fastq.gz]
[ctrl_1, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/ctrl_1.fastq.gz]
[ctrl_3, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/ctrl_3.fastq.gz]
[A375_1, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A375_1.fastq.gz]
[A375_2, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A375_2.fastq.gz]
[A375_3, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A375_3.fastq.gz]
[A549_1, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A549_1.fastq.gz]
[A549_2, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A549_2.fastq.gz]
[A549_3, /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/A549_3.fastq.gz]
```

Now we can use the format `tuple val(base), file(reads) from ch_reads` for kallisto alignment.

Read [here](https://www.nextflow.io/docs/latest/faq.html?highlight=simplename#how-do-i-get-a-unique-id-based-on-the-file-name) for descriptions on how to get unique IDs based on filenames using `simpleName` and `baseName`. 

## Nextflow script
Go to the main nextflow script available at this link: [https://github.com/BarryDigby/barrydigby.github.io/blob/master/RNA-Seq/main.nf](https://github.com/BarryDigby/barrydigby.github.io/blob/master/RNA-Seq/main.nf).

Save the contents as `main.nf` and invoke the help message by running `nextflow run main.nf --help`.  

> Note: the fastq files for the practical are located at: /data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/
