---
title: Walkthrough + Exercise
layout: page
permalink: /Variant_Calling/Nextflow
---

Due to time constraints, you are tasked with writing a nextflow script for the following processes:
1. Read Alignment
2. Mark Duplicates

I will provide a full nextflow script of the BASH workflow for you to run when you have completed the exercise.

***

# Reusable Channels
In the BASH workflow, you may have noticed nearly every step required the reference genome or the exome interval list file. In `nextflow` we typically intialise files with a channel, however usually channels can only be used once.

### Error example
Save the script as `non-reusable.nf` and run it -> `nextflow run non-reusable.nf`.

```bash
#!/usr/bin/env nextflow

Channel
	.fromPath("/data/MSc/2020/MA5112/Variant_Calling/reference/GRCh37.fasta")
	.set{ fasta_ch }

process A{
    	echo true

    	input:
    	file(fasta) from fasta_ch

    	script:
    	"""
    	printf "Process A, fasta file: $fasta\n"
    	"""
}

process B{
    	echo true

    	input:
    	file(fasta) from fasta_ch

    	script:
    	"""
    	printf "Process B, fasta file: $fasta\n"
    	"""
}
```

Note the error message:

> Channel `fasta_ch` has been used twice as an input by process `B` and process `A`

***

### `.into{}` example
One way to overcome this is to put the file into multiple channels using `.into{}` instead of `.set{}`.

Save the script as `into-reusable.nf` and run -> `nextflow run into-reusable.nf`

```bash
#!/usr/bin/env nextflow

Channel
	.fromPath("/data/MSc/2020/MA5112/Variant_Calling/reference/GRCh37.fasta")
	.into{ fasta_ch1; fasta_ch2 }

process A{
    	echo true

    	input:
    	file(fasta) from fasta_ch1

    	script:
    	"""
    	printf "Process A, fasta file: $fasta\n"
    	"""
}

process B{
    	echo true

    	input:
    	file(fasta) from fasta_ch2

    	script:
    	"""
    	printf "Process B, fasta file: $fasta\n"
    	"""
}
```

Process A and B can now run because the file has been put into two unique channels.

***

### `.value(file())` example
Placing a file into multiple channels is fine but is considered poor practice when writing large pipelines. For the variant calling workflow, we would have to place the fasta reference genome file into multiple channels, with each process hardcoded to a specific `fasta_ch(n)`. This is cumbersome to keep track of!

To create a reusable channel we are going to use `Channel.value(file())` to tell nextflow to stage the file path as a value and intialise it as a file. We can reuse this throughout the pipeline.

Save the below script as `value-reusable.nf` and run -> `nextflow run value-reusable.nf`

```bash
#!/usr/bin/env nextflow

Channel
	.value(file("/data/MSc/2020/MA5112/Variant_Calling/reference/GRCh37.fasta"))
	.set{ fasta_ch }

process A{
    	echo true

    	input:
    	file(fasta) from fasta_ch

    	script:
    	"""
    	printf "Process A, fasta file: $fasta\n"
    	"""
}

process B{
    	echo true

    	input:
    	file(fasta) from fasta_ch

    	script:
    	"""
      printf "Process B, fasta file: $fasta\n"
    	"""
}
```

We can see that the channel `fasta_ch` is reusable across processes.

***

# Staging files for Variant Calling
I have started you off by staging the input files required for the variant calling workflow. The Index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`, `.fai`, `.dict`) are all in a channel called `index_ch`. You will never need to specifically call one of these files, they just have to be present in the nextflow workdir for certain processes like read alignment (requires BWA indices) or BQSR (requires SAMtools `.fai` and Sequence Dictironary `.dict` files).

Other files such as `dbSNP`, `Mills_KG` or `exome_intervals` do need to be called specifically for flags. For this reason, they are in their own channel.

Save the script as `stage_inputs.nf` and run -> `nextflow run stage_inputs.nf`

```bash
#!/usr/bin/env nextflow

/*
 * Stage PE Sequencing Reads
 */

reads_ch = Channel.fromFilePairs("/data/MSc/2020/MA5112/Variant_Calling/reads/*_r{1,2}.fastq.gz")

/*
 * Stage Input files for Variant Calling:
 */

fasta_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/reference/GRCh37.fasta"))
exome_interval_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/assets/exome.bed.interval_list"))
mills_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/assets/Mills_KG_gold.indels.b37.vcf.gz"))
dbsnp_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/assets/dbsnp_138.b37.vcf.gz"))
index_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/reference/*GRCh37.{dict,fasta.}*"))

/*
 * I will show you how to call each of these files
 */

process A {
    	echo true

    	input:
    	tuple val(base), file(reads) from reads_ch
    	file(index) from index_ch
    	file(exome) from exome_interval_ch
    	file(dbsnp) from dbsnp_ch
    	file(mills) from mills_ch

    	script:
    	"""
    	echo $base, $reads
    	echo $index
    	echo $exome
    	echo $dbsnp
    	echo $mills
    	"""
}

/*
 * Check they can be reused.. (not the fastq pairs)
 */

index_ch.view()
exome_interval_ch.view()
dbsnp_ch.view()
mills_ch.view()
```


***

# Exercise
You are required to build on the `stage_inputs.nf` script and perform read alignment + markduplicates.

Refer to the full pipeline we ran for guidance. Below is a skeleton script to help:

```bash
#!/usr/bin/env nextflow

params.outdir = "./"

/*
 * Stage PE Sequencing Reads
 */

reads_ch = Channel.fromFilePairs("/data/MSc/2020/MA5112/Variant_Calling/reads/*_r{1,2}.fastq.gz")

/*
 * Stage Input files for Variant Calling:
 */

fasta_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/reference/GRCh37.fasta"))
exome_interval_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/assets/exome.bed.interval_list"))
mills_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/assets/Mills_KG_gold.indels.b37.vcf.gz"))
dbsnp_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/assets/dbsnp_138.b37.vcf.gz"))
index_ch = Channel.value(file("/data/MSc/2020/MA5112/Variant_Calling/reference/*GRCh37.{dict,fasta.}*"))

/*
 * ALIGN READS
 */

process bwa {
      publishDir "${params.outdir}/bwa_aln", mode:'copy'

      input:


      output:


      script:
      """

      """
}

/*
 * MARK DUPLICATES
 */

// Use the output BAM file as an input for this process

process markduplicates{
      publishDir "${params.outdir}/markdups", mode:'copy'

      input:


      output:


      script:
      """

      """
}
```

***

# Full script
If you want to see the full pipeline in nextflow, it is available (variant_calling.nf) at the following [repository](https://github.com/BarryDigby/barrydigby.github.io/tree/master/Variant_Calling). Note I used `.getval()` instead of `Channel.value(file())`, slightly outdated syntax but it still works.

Save the script to your own directory and run it by calling:

```bash
nextflow -bg -q run variant_calling.nf \
--outDir $(pwd) \
--analysisDir "/data/MSc/2020/MA5112/Variant_Calling" \
-with-singularity /data/MSc/2020/MA5112/Variant_Calling/contaier/germline_vc.img
```
