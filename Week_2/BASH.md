---
title:
layout: page
permalink: /Week_2/BASH
---

# Exome Sequencing Analysis
We are going to analyse a subset of a patients whole exome sequencing data to identify variants in protein coding genes, attempting to elucidate possible genetic diseases.

# Workflow

### Genome Index
As with any analysis, the reference genome must be indexed to quickly extract alignments overlapping particular genomic regions.

```bash
bwa index -a bwtsw GRCh38.fa
```

* `-a`: Construction algorithm (bwtsw, is or rb2). For large genomes, specify `-a bwtsw`. If you are unsure, omit this flag and `bwa` will determine the correct algorithm to use.

##### Output
Indexing using `bwa` produces 5 files: `*.amb`, `*.ann`, `*btw`, `*.pac` & `*.sa`.

```bash
samtools faidx GRCh38.fa
```

##### Output
Samtools generates a `*.fai` file (<strong>fa</strong>sta <strongi</strong>ndexed). 
