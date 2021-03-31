---
title: Kallisto Pipeline
layout: page
permalink: /RNA-Seq/pipeline
---

A description of the `kallisto` pipeline has been provided for you below.

There is no need to run this for this weeks tutorial, it has been posted as a learning resource.

# 1. Indexing
Kallisto requires an indexed transcriptome file for downstream quantification.

##### *Inputs*
- Reference [transcriptome file](http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz).

```bash
kallisto index -i GRCh38.idx Homo_sapiens.GRCh38.cdna.all.fa
```

##### *Flags*
* `-i`: Output filename

##### *Outputs*
Kallisto index produces an indexed transcriptome file

***

# 2. Quantification
Kallisto can perform quantification using either single-end or paired-end fastq files.

##### *Inputs*
- Indexed transcriptome
- FASTQ files

##### Paired-end
```bash
kallisto quant \
         -i Genome Index file \
         -t 2 \
         -o outDir/ \
         --bias \
         FASTQ_1, FASTQ2
```

##### Single-end
```bash
kallisto quant \
         --single \
         -l 200 \
         -s 30 \
         -i Genome Index file \
         -t 2 \
         -o outDir/ \
         --bias \
         FASTQ
```

##### *Flags*
* `-i` Indexed genome file from `kallisto index`
* `--single` Indicate input is single-end reads (requires `-l` and `-s`)
* `-t` *n* threads to use
* `-o` Output directory of the sample, containing 3 files. Do not name the directory `outDir`, name it according to the sample name e.g CTRL_2.fastq should have the directory name CTRL_2/ etc.
* `-l` Estimated average fragment length
* `-s` Estimated standard deviation of fragment length
* `--bias` Perform sequence based bias correction

> Note: Estimated fragment lengths and standard deviation must be retrieved from the sequencing center. If using public data, try to find this information online. Using incorrect values will greatly influence the outputs. 200 and 30 are typically used for average fragment length and standard deviation, respectively. 

##### *Outputs*
Kallisto `quant` will output a directory for each sample containing:

- `abundance.h5`
- `abundance.tsv`
- `run_info.json`

The two abundance files contain transcript quantification information. For downstream analysis we can use either the `.h5` file or the `.tsv` file. The difference between the `.h5` and `.tsv` file is the `.h5` file contains bootstrapping information if the option was specified during kallisto quant. This is for downstream analysis using Sleuth (not covered in this tutorial).

We did not specify the `-bootstrap option`, thus we can use either the `.h5` or `.tsv` file for analysis in R as they are identical.
