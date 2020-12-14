---
title: Anaconda
layout: page
permalink: /Week_1/Conda
---

Topics covered:
- [Conda Install] (#install)
- [Conda yaml file] (#yaml)


Anaconda is an open source distribution of packages for Linux, Mac OS and Windows. Tools are installed using `conda` package and environment manager, creating a clean slate environment in which it attempts to install packages without conflicts. Conda offers a one stop shop for genomic tools and packages, which can be installed using a simple `conda install` command.

## Conda Install {#install}

***

*Please run this locally on your own laptop*

***

Let's first activate the (base) environment for conda:

```bash
conda activate base
```

You should notice the `(base)` prefix before your username on the terminal. Now install the package `fastqc`. Append the channel in which to search for the package. Looking at the Anaconda repository, we can see that it belongs to `bioconda` https://anaconda.org/bioconda/fastqc. If we want to pin the tool version number, we can install it using the command `conda install bioconda::fastqc=0.11.9`

```bash
conda install bioconda::fastqc

or

conda install -c bioconda fastqc
```

 Test the package works by prompting the help documentation with the `-h` flag:

 ```bash
 fastqc -h
 ```

***

## Conda yaml file (#yaml)

We have seen how simple it is to install tools using `conda` by using the `conda install command`. However in reality we will want to install multiple tools at once for an analysis and create a clean environment for the tools. This can be simplified using a `.yml` file.

As this tutorial covers quality control in sequencing reads, we will need a tool to trim and remove adapters in addition to the fastqc tools. Below is an example of a `.yml` file to use for creating a quality control environment:

```
name: QC
channels:
  - agbiome
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - fastqc
  - multiqc
  - bbtools
```
