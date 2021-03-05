---
title: Nextflow
layout: page
permalink: /RNA-Seq/Nextflow
---

# Walkthrough + Exercise
Before beginning the exercise I will demonstrate using the `map` and `combine` operators in nextflow.

Save the below script as `test1.nf` and run `nextflow run test1.nf`

```nextflow
#!/usr/bin/env nextflow

// Initialise reference TX file channel
transcriptome_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reference/*.fa")

// Initialise reads channel (reads are single end)
reads_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reads/*.fastq.gz")

transcriptome_ch.view()

reads_ch.view()
```

We can see each channel holds the path to each file.

***

Due to the fact out reads are single-end, we used `fromPath()` and not `fromFilePairs()`. We are at a slight disadvantage here, as the files are not automatically placed in a named tuple.

We will use the `map` operator to assign the files `simpleName` which strips all text after the first period character. (This is in contrast to `baseName`, which will remove all text after the **last** period character: `ctrl_1.fastq.gz` becomes `ctrl_1.fastq`).

Save the below script as `test2.nf` and run `nextflow run test2.nf`

```nextflow
#!/usr/bin/env nextflow

// Initialise reference TX file channel
transcriptome_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reference/*.fa")

// Initialise reads channel (reads are single end)
reads_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reads/*.fastq.gz")

transcriptome_ch.view()

// use the map operator
reads = reads_ch.map{ file -> [file.simpleName, file]}

reads.view()
```

Notice that all fastq reads have been placed in a tuple with their simpleName. This means we can use `tuple val(base) file(read) from reads`, which is useful for dynamically assigning the output names of files.

***

Lets use the two chanels in a process. We are going to `echo` the files, however this could be a `kallisto` command. The idea here is to make sure all files are being used by the process properly!

Save the below script as `test3.nf` and run `nextflow run test3.nf`

```nextflow
#!/usr/bin/env nextflow

// Initialise reference TX file channel
transcriptome_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reference/*.fa")

// Initialise reads channel (reads are single end)
reads_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reads/*.fastq.gz")

// use the map operator, assign to new channel
reads = reads_ch.map{ file -> [file.simpleName, file]}

process foo {

	echo true

	input:
	tuple val(base), file(read) from reads
	file(idx) from transcriptome_ch

	output:
	stdout to out

	script:
	"""
	echo "File Prefix: $base, File: $read"
	echo ""
	echo "Transcriptome file: $idx"
	"""
}
```

**Notice that only one fastq file is being used**

We wanted the process to submit a job for each fastq file and run in parallel - one of the most attractive reasons to use nextflow in the first place.

***

To remedy the above problem, we will use the `combine` operator. This will add the transcriptome index file to each read tuple: [ base, [ fastq, index]]. When calling them in `inputs:`, the channel will acknowledge that the index file came from a different channel and can thus be called with its own `file(idx)` call. 

Save the below script as `test4.nf` and run `nextflow run test4.nf`

```nextflow
#!/usr/bin/env nextflow

// Initialise reference TX file channel
transcriptome_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reference/*.fa")

// Initialise reads channel (reads are single end)
reads_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reads/*.fastq.gz")

// use the map operator, assign to new channel
reads = reads_ch.map{ file -> [file.simpleName, file]}

process foo {

	echo true

	input:
	tuple val(base), file(read), file(idx) from reads.combine(transcriptome_ch)

	output:
	stdout to out

	script:
	"""
	echo "File Prefix: $base, File: $read"
	echo ""
	echo "Transcriptome file: $idx"
	"""
}
```

# Exercise
Refer to the `pipeline` section of this weeks tutorial for tips on how to run kallisto. The previous tips using `map` and `combine` should provide you with eveything you need to execute the analysis.

***

Below is the skeleton script. Complete the script and ask us to check it before running it in your own directory.

```
nextflow -bg -q run kallisto.nf \
--outDir $(pwd) \
-with-singularity /data/MSc/2020/MA5112/RNA-Seq/container/RNA-Seq.img
```


```nextflow
#!/usr/bin/env nextflow

// Initialise reference TX file channel
transcriptome_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reference/*.fa")

// Initialise reads channel (reads are single end)
reads_ch = Channel.fromPath("/data/MSc/2020/MA5112/RNA-Seq/reads/*.fastq.gz")

// use the map operator, assign to new channel
reads = reads_ch.map{ file -> [file.simpleName, file]}

process Index{

    publishDir "$params.outDir/index/", mode:'copy'

    input:


    output:


    script:
    """

    """
}

process Quantification{

    publishDir "$params.outDir/quant/", mode:'copy'

    input:


    output:


    script:
    """

    """
}
```
