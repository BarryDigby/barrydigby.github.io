---
title:
layout: page
permalink: /Week_2/BASH
---

# Exome Sequencing Analysis
We are going to analyse a subset of a patients whole exome sequencing data to identify variants in protein coding genes, attempting to elucidate possible genetic diseases.

# Workflow

### 1.Genome Index
As with any analysis, the reference genome must be indexed to quickly extract alignments overlapping particular genomic regions.

### BWA Index
```bash
bwa index -a bwtsw GRCh38.fa
```

* `-a`: Construction algorithm (bwtsw, is or rb2). For large genomes, specify `-a bwtsw`. If you are unsure, omit this flag and `bwa` will determine the correct algorithm to use.

##### Output
Indexing using `bwa` produces 5 files: `*.amb`, `*.ann`, `*btw`, `*.pac` & `*.sa`.

***

### Samtools Index
```bash
samtools faidx GRCh38.fa
```
##### Output
Samtools generates a `*.fai` file (<strong>fa</strong>sta <strong>i</strong>ndexed).

***

### 2. Map Reads
Map the reads to the reference genome using `bwa mem`, and pipe the output to `samtools sort`.

```bash
bwa mem \
    -K 10000000 \
    -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' \
    -t 8 \
    GRCh38.fa \
    ../reads/subsample_r1.fastq.gz \
    ../reads/subsample_r2.fastq.gz \
    | samtools sort - > subsample.bam
```

* `-K`: process INT input bases in each batch regardless of nThreads (for reproducibility)
* `-R`: Read Group identifier. This information is key for downstream GATK functionality, GATK will not work without a read group tag.
* `-t`: Number of threads

The script above passes the output of `bwa mem` (`subsample.sam`) directly to `samtools sort` by using the pipe command (`|`). Samtools uses a dash `-` to indicate the incoming file from the previous process. When working with full sized samples, allocate more threads to samtools like so `samtools sort --threads 8 - > subsample.bam`.
