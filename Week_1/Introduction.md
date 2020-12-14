---
title: Week 1 Introduction
layout: page
permalink: /Week_1/Introduction
---

This tutorial covers creating a container for quality control of reads and using the container in a nextflow script to perform QC of downloaded reads.

## Container Creation
Before constructing a container, you need to know exactly which tools are going to be used for the analysis. For the quality control of sequencing reads, we will need `FastQC` and `MultiQC` to generate `.html` summary statistics of the sequencing reads. We will also use `bbtools` to perform adapter trimming and low quality read removal.

When creating a Docker or Singularity container, one can specify the tools to download in the Dockerfile or Singularity definition file. This requires one to install each tool manually and append it to the `$PATH` of the container.

Alternatively, a much easier way to download tools and append them to the `$PATH` is to specify a `.yml` file containing the tools to download via Anaconda repository. We will create a conda environment within the conatiner and append its executables to the containers `$PATH`. Creating containers in this way is safer, as conda will warn you if tools have conflicting dependencies. It is also highly reproducible, all that are needed to create the container is a `.yml` file and a `Dockerfile` file.

##### Conda .yml file
Below is an example of a conda `.yml` file used to create an environment called `week1`

```
name: week1
channels:
  - agbiome
  - bioconda
  - conda-forge
dependencies:
  - bbtools
  - fastqc
  - multiqc
```

To create a conda environment with the above tools installed, run the following code:

```bash
conda env create -f week1.yml
```

A conda environment has now been created for you. Lets activate the environment and check the tools were correctly installed by calling their help messages:

```bash
conda activate week1
fastqc -h
multiqc -h
bbduk.sh -h
```

To help explain how Docker uses the conda environment, we will first check where the conda environment we just created is installed. While `week1` environment is active, type `whereis bbduk.sh`.

```bash
whereis bbduk.sh
bbduk: /home/barry/anaconda3/envs/week1/bin/bbduk.sh
```

We can see the path where the executable `bbduk.sh` is installed. When conda activates an environment, it appends the environment bin (/home/barry/anaconda3/envs/week1/bin/) to your `$PATH`, so you can use tools in the environment without having to specify a full `$PATH`.

Docker takes advantage of this by creating a conda environment using the `.yml` file specified, and appends the environments path to the containers path. This will be shown in the next section.

##### Dockerfile
To create a Docker container, a `Dockerfile` is required in the directory. This is a file that specifies the build and contents of the container. Make sure the `week1.yml` file specified above is in your current working directory. Add the `Dockerfile` contents given below:

```bash
FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="MA5112 week1 tutorial"

WORKDIR ./
COPY week1.yml ./
RUN conda env create -f week1.yml && conda clean -a
ENV PATH /opt/conda/envs/week1/bin:$PATH
```

The Docker commands used in the `Dockerfile`:
- `FROM`: Creates a layer for the image. We are using `nextflow` `nfcore` base as it comes pre-installed with conda.
- `LABEL`: Description of the image.
- `WORKDIR`: Change the work directory of the container. We are pointing to the current directory `./` with the `week1.yml` file.
- `COPY`: Copy contents from the `WORKDIR` to the container. We want to use the `week1.yml` file for conda, so we copy it locally to the container.
- `RUN`: Execute a command.
- `ENV`: Edit the containers `$PATH` environment variable in its `~/.bashrc`.
