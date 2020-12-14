---
title: Overview
layout: page
permalink: /RNA-Seq/Overview
---

Before constructing a container and nextflow script for the analysis I will walk through what is typically involved in an RNA-Seq analysis. This will be done using the command line which we will end up utilising in the nextflow script later.

- [Required Files](#required)
- [Dataset](#data)
- [Quality Control](#quality)
- [Deleting files](#deleting)

## Required Files {#required}
RNA-Seq comes in two flavors - gene quantification and transcriptome assembly. Gene expression quantification is the more common of the two - we are interested in the levels of gene expression by way of counts. Transcript assembly can be either reference free (for species that are not well annotated) or reference guided to provide insights into novel transcripts that are a by product of gene fusions or non-canonical splicing.

In the tutorial we will be covering gene expression quantification. The reference file required for this contains transcript sequences, as we are only interested in cDNA and do not require the full genome sequence. Personally I prefer using `gencode`, however the same file is available at `ENSEMBL`.

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz && gunzip gencode.v36.transcripts.fa.gz && mv gencode.v36.transcripts.fa GRCh38.cDNA.fa
```

Now we have a transcript reference file, `GRCh38.cDNA.fa` which will be the reference file for the analysis.

It is available at the path `/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Reference/GRCh38.cDNA.fa`

***

## Dataset {#data}
<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/rna-seq/rnaseq-data.png" width="100%" height="100%"/>
</center>

The researcher has captured exosomes from human malignant melanoma cell line (A375) and the human lung cancer cell line (A549) and added them to primary human epidermal melanocytes, (NHEM-c cells). There are triplicates of each cell line:

- NHEM-c (Control)
- NHEM-c + A549 exosomes (Lung)
- NHEM-c + A375 exosomes (Melanoma)

The researcher is interested in characterising the transcriptional profile of each transfected cell line, in an effort to elucidate the tumor pathways mediated by A375 & A549 tumor exosomes.

## Quality Control {#quality}
As with all NGS assays, the quality of the reads must be assessed and checked for adapter contamination.
