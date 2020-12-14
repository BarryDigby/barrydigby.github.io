---
title: Docker
layout: page
permalink: /Week_1/Docker
---

***

*Please create an account on Dockerhub before continuing!*

[https://hub.docker.com/](https://hub.docker.com/)

***

Jump to:
 - [Dockerfile](#dockerfile)
 - [Create Container](#create)
 - [Push to Dockerhub](#push)

We are going to take advantage of the conda `.yml` file created in the previous step to create a Docker container with the environment `QC` appended to it's `$PATH`.

***

# Dockerfile {#dockerfile}
To create a Docker container, we need to write a `Dockerfile` containing the instructions on which layer and packages to use. Below is an example of a `Dockerfile` using the conda `.yml` file created in the the Conda section of this weeks tutorial:

```bash
FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="MA5112 week1 tutorial"

WORKDIR ./
COPY week1.yml ./
RUN conda env create -f week1.yml && conda clean -a
ENV PATH /opt/conda/envs/week1/bin:$PATH
```

The commands used in the `Dockerfile`:
- `FROM`: Creates a layer for the image. We are using `nextflow` `nfcore` base as it comes pre-installed with conda.
- `LABEL`: Description of the image.
- `WORKDIR`: Change the work directory of the container. We are pointing to the current directory `./` containing the `week1.yml` file.
- `COPY`: Copy contents from the `WORKDIR` to the container. We want to use the `week1.yml` file for conda, so we copy it to the container.
- `RUN`: Execute a command.
- `ENV`: Edit the containers `$PATH` environment variable in its `~/.bashrc`.

## Create Container {#create}
*Personal Dockerhub required* -- **Please do not push to my dockerhub account**.

***

Make sure the `Dockerfile` and `.yml` file are in your current working directory.

Go to dockerhub, login and create a new repository "week1". We will tag the container using your username and repository:

```bash
docker build -t barryd237/week1:test .
```


## Push to Dockerhub {#push}
The container has been created locally, now let's push it to dockerhub:

```bash
docker run barryd237/week1:test
docker push barryd237/week1:test
```

The container is now pushed to my dockerhub account at [https://hub.docker.com/repository/docker/barryd237/week1](https://hub.docker.com/repository/docker/barryd237/week1).

***

## Try it yourself

Try publishing the container to your own dockerhub profile, following the structure of `docker build -t DOCKERUSERNAME/REPO:TAG`.

The `tag` can be anything, however in a situation where you are working on developing a pipeline use tags such as `:dev` to delineate between development versions and published versions.
