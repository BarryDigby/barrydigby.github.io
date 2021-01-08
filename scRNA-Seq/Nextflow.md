---
title: Nextflow
layout: page
permalink: /scRNA-Seq/Nextflow
---

# Walkthrough
The analysis using `kallisto` and `bustools` uses directory as outputs and inputs to downstream processes which we have not encountered before.

I will demonstrate how these work with some proof of concept scripts below.

***

Save the below script as `test.nf` and run `nextflow run test.nf`

```bash
#!/usr/bin/env nextflow


// create a directory in a process
// add some empty files

process foo{

	output:
	file("folder") into foo_out

	script:
	"""
	mkdir -p folder/
	touch file{1..3}.txt
	mv *.txt folder/
	"""
}

// place foo_out into two new channels
(foo1, foo2) = foo_out.into(2)

//view the channel
foo1.view()

process bar{

	echo true

	input:
	file(dir) from foo2

	output:
	file(folder) into bar_out

	script:
	file1 = "${dir}/file1.txt"
	"""
	echo "$file1"
	"""
}

process baz{

	echo true

	input:
	file(dir) from bar_out

	output:
	stdout to out

	script:
	"""
	for i in $dir/*; do

		printf "\$i\n"

	done
	"""
}
```

***

# Analysis
There will be no exercise for the scRNA-Seq workflow as it is slightly complex.

Save the below script as `kallisto_bustools.nf` and run in your directory:

```bash
nextflow -bg -q \
run kallisto_bustools.nf \
--outDir $(pwd) \
-with-singularity /data/MSc/2020/MA5112/scRNA-Seq/container/scRNA.img
```

```nextflow
#!/usr/bin/env nextflow

// FASTQ reads
params.reads = "/data/MSc/2019/1k_pbmc_protein_v3_gex_fastqs/*_R{1,2}_*"
Channel
	.fromFilePairs(params.reads)
	.set{reads_ch}

// 10x Genomics whitelist
params.whitelist = Channel.fromPath("/data/MSc/2020/MA5112/scRNA-Seq/assets/10xv3_whitelist.txt").getVal()

// Transcripts to Gene
params.tx2gene = Channel.fromPath("/data/MSc/2020/MA5112/scRNA-Seq/assets/tx2gene.txt").getVal()

// Kallisto Index file (pre-made to save time)
params.idx = Channel.fromPath("/data/MSc/2020/MA5112/scRNA-Seq/reference/Homo_sapiens.cDNA.idx").getVal()

// Experiment chemistry (v2 or v3)
params.chemistry = "10xv3"

/*
  Kallisto Bustools
*/

process kallisto{

	publishDir "$params.outDir/analysis/raw", mode:'copy'

	input:
	tuple val(base), file(reads) from reads_ch
	file(idx) from params.idx

	output:
	file("bus_output") into kallisto_bus_sort
	file("kallisto.log") into log_out

	script:
	"""
	kallisto bus \
	-i $idx \
	-o bus_output/ \
	-x ${params.chemistry} \
	-t 2 \
	$reads | tee kallisto.log
	"""
}


process bustools_correct_sort{

	publishDir "$params.outDir/analysis/sorted", mode:'copy'

	input:
	file(bus) from kallisto_bus_sort
	file(whitelist) from params.whitelist

	output:
	file(bus) into kallisto_corrected

	script:
	correct = "bustools correct -w $whitelist -o ${bus}/output.corrected.bus ${bus}/output.bus"
	sort_file = "${bus}/output.corrected.bus"
	"""
	$correct
	mkdir -p tmp
	bustools sort \
	-T tmp/ \
	-t 2 \
	-o ${bus}/output.sort.bus \
	$sort_file
	"""
}



process bustools_count{

	publishDir "$params.outDir/analysis/counts", mode:'copy'

	input:
	file(bus) from kallisto_corrected
	file(tx2gene) from params.tx2gene

	output:
	file("output.sort.txt") into out

	script:
	"""
	bustools text \
	-o output.sort.txt \
	${bus}/output.sort.bus
	"""
}
```
