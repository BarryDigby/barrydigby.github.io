---
title: Week 1 Introduction
layout: page
permalink: /Week_1/Introduction
---

The week 1 tutorial will show you how to use `Anaconda` to construct a `Docker` container in a reproducible manner, create an image of the container using `Singularity` and use the image in a `nextflow` script to perform qualtiy control on sequencing reads.

Topics covered this week:

- [1. Conda](http://barrydigby.github.io/Week_1/Conda)
- [2. Docker](http://barrydigby.github.io/Week_1/Docker)
- [3. Singularity](http://barrydigby.github.io/Week_1/Singularity)


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


## Docker Hub
Please create you own Docker hub account before proceeding with the tutorial. We will publish our created containers to dockerhub which will enable automatic linting and allows us to pull the image from the containers docker page.

## Creating the container
Now that we have a docker account, `Dockerfile` and `.yml` file we can begin to create the container.

```bash
docker build -t barryd237/week1:test
```

The above command creates the docker container, tagging (`-t`) it under the account barryd237, in the week1 repository with the label test.

## Publishing to Docker
Now that the container has been creted, lets push it to docker hub.

```bash
docker run barryd237/week1:test
docker push barryd237/week1:test
```
 Check the container has been pushed to your dockerhub profile. My container is available at: https://hub.docker.com/r/barryd237/week1.

## Create an image
The docker container has been created and pushed to dockerhub. Next we will use singularity to pull the docker container and write it to an image file.

```bash
singularity pull --name week1.img docker://barryd237/week1:test
```

You will have an image in your directory called `week1.img`. We can shell into the container using `singularity shell -B $(pwd) week1.img`
.
