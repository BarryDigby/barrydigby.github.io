---
title: Kallisto Pipeline
layout: page
permalink: /RNA-Seq/pipeline
---

A description of the `kallisto` pipeline has been provided for you below.

There is no need to run this for this weeks tutorial, it has been posted as a learning resource.

# 1. Indexing
Kallisto requires an indexed genome file for downstream quantification.

##### *Inputs*
- Reference Genome

```bash
kallisto index -i GRCh37.idx GRCh37.fa
```

##### *Flags*
* `-i`: Filename for the kallisto index to be constructed

##### *Outputs*
Kallisto index produces an indexed genome file with the filename specified by `-i`

***

# 2. Quantification
Kallisto can perform quantification using either single-end or paired-end fastq files.

##### *Inputs*
- Indexed Genome
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
* `-o` Output directory of the sample, containing 3 files
* `-l` Estimated average fragment length
* `-s` Estimated standard deviation of fragment length
* `--bias` Perform sequence based bias correction

##### *Outputs*
Kallisto `quant` will output a directory for each sample containing:

- `abundance.h5`
- `abundance.tsv`
- `run_info.json`

The two abundance files contain transcript quantification information. For downstream analysis we can use either the `.h5` file or the `.tsv` file. The difference between the `.h5` and `.tsv` file is the `.h5` file contains bootstrapping information if the option was specified during kallisto quant. This is for downstream analysis using Sleuth (not covered in this tutorial).

We did not specify the `-bootstrap option`, thus we can use either the `.h5` or `.tsv` file for analysis in R as they are identical.
