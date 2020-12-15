---
title:
layout: page
permalink: /Week_1/Docker
---

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/week1/docker.png" width="100%" height="100%"/>
</center>

> "A container is a standardized unit of software. It helps to distribute the software regardless of the setup. It packages up code and all its dependencies so the application runs quickly and reliably from one computing environment to another."


Jump to:
 - [Dockerfile](#dockerfile)
 - [Create Image](#create)
 - [Push to Dockerhub](#push)

# Dockerfile {#dockerfile}
To create a Docker image, we need to write a `Dockerfile` containing the instructions on which layer and packages to use.

Below is an example of a Dockerfile using the `week1.yml` file created in the the Anaconda section of this weeks tutorial:

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
- `FROM`: Creates a layer for the image. We are using `nfcore` base, a Ubuntu distro with conda pre-installed.
- `LABEL`: Description of the image.
- `WORKDIR`: Specify the local working directory for the image.
- `COPY`: Copy contents from the `WORKDIR` to the image. We want to use the `week1.yml` file for conda, so we copy it to the image.
- `RUN`: Execute a command.
- `ENV`: Append `$PATH` environment variables in images `~/.bashrc`.

# Create Image {#create}
> "A Docker image is an immutable file that contains the source code, libraries, dependencies, tools, and other files needed for an application to run."

We will build an image and tag it using your dockerhub details.

### To Do:
1. Create Dockerhub repository 'week1'
2. Make sure `Dockerfile` & `week1.yml` are in current directory.
3. Build Docker image
4. List image using `docker images`.

```bash
docker build -t USERNAME/REPO:TAG .
```

# Push to Dockerhub {#push}
The image has been created locally, now let's push it to dockerhub:

```bash
docker run USERNAME/REPO:TAG
docker push USERNAME/REPO:TAG
```

Go to your Dockerhub account and check the image has been stored under the repository week1.

*N.B: The `TAG` can be anything, however in a situation where you are developing a pipeline it is standard practice to use the dev tag to test new features before rolling them out as a new version.*
