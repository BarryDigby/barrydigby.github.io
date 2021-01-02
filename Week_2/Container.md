---
title: Variant Calling Container
layout: page
permalink: /Week_2/Container
---

Below are the resources required to construct a variant calling container. I have included several extra variant calling tools in the container (`Strelka`, `Mantra`, `Google DeepVariant`, `GATK`, `FreeBayes`) as some of you might want to try different variant calling tools for your thesis project.

### Jump to
***
- [yml](#yml)
- [Dockerfile](#Dockerfile)
- [Singularity](#Singularity)

# yml {#yml}
```
name: Germline_VC
channels:
  - bioconda
  - conda-forge
  - defaults
  - agbiome
dependencies:
  - agbiome::bbtools
  - conda-forge::markdown=3.1.1
  - conda-forge::pymdown-extensions=6.0
  - conda-forge::pygments=2.5.2
  - bioconda::ascat=2.5.2
  - bioconda::bcftools=1.9
  - bioconda::bwa=0.7.17
  - bioconda::cancerit-allelecount=4.0.2
  - bioconda::cnvkit=0.9.6
  - bioconda::control-freec=11.5
  - bioconda::ensembl-vep=99.2
  - bioconda::fastqc=0.11.9
  - bioconda::freebayes=1.3.2
  - bioconda::gatk4-spark=4.1.7.0
  - bioconda::genesplicer=1.0
  - bioconda::htslib=1.9
  - bioconda::manta=1.6.0
  - bioconda::msisensor=0.5
  - bioconda::multiqc=1.8
  - bioconda::qualimap=2.2.2d
  - bioconda::samtools=1.9
  - bioconda::snpeff=4.3.1t
  - bioconda::strelka=2.9.10
  - bioconda::tiddit=2.7.1
  - bioconda::trim-galore=0.6.5
  - bioconda::vcfanno=0.3.2
  - bioconda::vcftools=0.1.16
  - conda-forge::pigz=2.3.4
  - conda-forge::r-ggplot2=3.3.0
  - conda-forge::pip=10.0.1
  - bioconda::deepvariant=0.10.0
  - bioconda::picard=2.18.7
  - conda-forge::lbzip2=2.5
  - conda-forge::google-cloud-sdk<243.0.0
  ```

# Dockerfile {#Dockerfile}
```
FROM nfcore/base:1.9
LABEL authors="Barry Digby" \
      description="Docker image containing all software requirements for Germline Variant Calling Analysis"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml python=2.7.15 && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/Germline_VC/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name Germline_VC > Germline_VC.yml
```

# Singularity {#Singularity}
You can create the image yourself and push to your own dockerhub account or you can download it from my dockerhub account.

```
singularity pull --name germline_vc.img docker://barryd237/germline_vc:latest
```
