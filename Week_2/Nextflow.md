---
title: Week 2 Exercise
layout: page
permalink: /Week_2/Nextflow
---

Due to time constraints, you are tasked with writing a nextflow script for the following processes:
1. Read Alignment
2. Mark Duplicates

***

## getVal()
In the BASH workflow, you may have noticed nearly every step required the reference genome. In `nextflow` we typically intialise files with a channel, however, as channels can only be used once we would have to create multiple channels for the reference genome file (something like `fasta_a`, `fasta_b`. `fasta_c` etc.) for each step in the workflow.

Instead we are going to use `getVal()` when specifying the parameter file path. When using the parameter in a process we have two choices:
1. Call the file on its own: `file(var_name) from params.var`
2. Combine files in a tuple: `tuple file(var1), file(var2) from Channel.value([params.var1, params.var2])`

Run the below nextflow script to gain some intuition: `nextflow run test.nf`

```nextflow
#!/usr/bin/env nextflow

// Reference Genome
params.fasta = Channel.fromPath("/data/MSc/2020/MA5112/Variant_Calling/reference/*fasta").getVal()
params.fai = Channel.fromPath("/data/MSc/2020/MA5112/Variant_Calling/reference/*fai").getVal()

// Process foo uses the reference genome only

process foo{
    	echo true

    	input:
    	file(fasta) from params.fasta

    	output:
   	stdout to foo_out

    	script:
    	"""
   	echo "process foo output"
    	echo "Reference Genome file:  $fasta"
    	"""
}


// process bar uses both the reference genome and the samtools index file
// Call them both in the same line using Channel.value() placing them in a tuple

process bar{
    	echo true

    	input:
    	tuple file(fasta), file(fai) from Channel.value([params.fasta, params.fai])

    	output:
    	stdout to bar_out

    	script:
    	"""
    	echo "process bar output"
    	echo "Reference genome file: $fasta"
    	echo "Samtools index file: $fai"
    	"""
}

// process baz modifies the index file (trivial example)
// to demonstrate examples of basic outputs

process baz{

	input:
	file(fai) from params.fai

	output:
	file("*.fai") into baz_output

	script:
	"""
	mv $fai new_index.fai
	"""
}

// In a nextflow process, you can only capture what the process
// created to an output channel (unless you use special flags)
// this is why GRCh37.fai is not captured by *.fai

// lets view the contents of the channel to confirm this:
baz_output.view()

/*
  FASTQ Example
*/

params.reads = "/data/MSc/2020/MA5112/Variant_Calling/reads/*_r{1,2}.fastq.gz"
Channel
 .fromFilePairs(params.reads)
 .into{reads_ch; reads_ch_print}

// set{a} for one channel
// into{a;b;c;...} for multiple channels
// view the tuple in channel
reads_ch_print.view()

process qux{
	echo true

	input:
	tuple val(base), file(reads) from reads_ch

	output:
	tuple val(base), file("*.fq.gz") into qux_out

	script:
	"""
	echo "tuple val(base) uses the common string in the fastq read pair name, dictated by the glob pattern in params.reads"
	echo "in our case, the variable base is: $base"
	echo ""
	echo 'The variable (dollar)reads stores r1 & r2 in the order they appear in the tuple'
	echo "$reads"

	# to access individual read pairs, use 0 based indexing
	# use ${base} to name outputs, appending the correct extension
	mv ${reads[0]} ${base}_r1.fq.gz
	mv ${reads[1]} ${base}_r2.fq.gz
	"""
}

qux_out.view()
```

***

# Exercise
I have initialised all of the parameters required to run the first 2 processes. Look back at the BASH workflow to see what files are required for each step.

*N.B* pass `--analysisDir` via the command line (`/data/MSc/2020/MA5112/Variant_Calling`) & `--outDir` when running the script. Assuming you are running the script in your personal directory on lugh (`/data/MSc/2020/username`), set this to `.`

```bash
nextflow run germline_vc.nf \
--analysisDir /data/MSc/2020/MA5112/Variant_Calling \
--outDir /data/MSc/2020/username \
-with-singularity /data/MSc/2020/MA5112/Variant_Calling/container/germline_vc.img
```

You must fill in all inputs/outputs/script body for both processes or the script will not run. When you think you have completed the script, let me know and I will check it before you run it.

```nextflow
#!/usr/bin/env nextflow

/*
Sequencing Reads Channel
*/

params.reads = "/data/MSc/2020/MA5112/Variant_Calling/reads/*_r{1,2}.fastq.gz"
Channel
 .fromFilePairs(params.reads)
 .set{reads_ch}

/*
Reference Genome Files
*/

params.fasta = Channel.fromPath("$params.analysisDir/reference/*fasta").getVal()
params.fai = Channel.fromPath("$params.analysisDir/reference/*.fasta.fai").getVal()
params.dict = Channel.fromPath("$params.analysisDir/reference/*.dict").getVal()

params.amb = Channel.fromPath("$params.analysisDir/reference/*fasta.amb").getVal()
params.ann = Channel.fromPath("$params.analysisDir/reference/*fasta.ann").getVal()
params.bwt = Channel.fromPath("$params.analysisDir/reference/*fasta.bwt").getVal()
params.pac = Channel.fromPath("$params.analysisDir/reference/*fasta.pac").getVal()
params.sa = Channel.fromPath("$params.analysisDir/reference/*fasta.sa").getVal()

/*
 Directories
*/

params.outDir = ""
params.analysisDir = ""

/*
================================================================================
                                 ALIGNMENT
================================================================================
*/

process MapReads{

publishDir path: "$params.outDir/analysis/bwa_aln", mode: "copy"

    input:


    output:


    script:
    readGroup = "@RG\\tID:sample_1\\tLB:sample_1\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:sample_1"
    """
    bwa mem \

    -R \"${readGroup}\" \




    """
}


/*
================================================================================
                                PROCESSING
================================================================================
*/


process MarkDuplicates{

    publishDir path: "$params.outDir/analysis/mark_dups", mode: "copy"

    input:


    output:


    script:
    """







    """
}
```
