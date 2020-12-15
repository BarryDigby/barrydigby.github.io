---
title:
layout: page
permalink: /Week_1/Conda
---

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/week1/anaconda_horizontal.png" width="100%" height="100%"/>
</center>

Anaconda is an open source distribution of packages for POSIX compute environments. Conda quickly installs, runs and updates packages and their dependencies. Conda easily creates, saves, loads and switches between environments on your computer, resolving dependency conflicts that might arise due to package requirements.

Jump to:
- [Install Packages](#install)
- [YAML file](#yaml)

# Install Packages {#install}
***
Let's first activate the base environment for conda:

```bash
conda activate base
```

You should see the `(base)` prefix before your username on the terminal.

Now let's install the package `fastqc`. First, look up the package in the Anaconda repository: [https://anaconda.org/bioconda/fastqc](https://anaconda.org/bioconda/fastqc).

```bash
conda install -c bioconda fastqc
```

*or*

```bash
conda install bioconda::fastqc
```

If we wanted to specify the version of the tool, we can 'pin' the version in the install command:

```bash
conda install bioconda::fastqc=0.11.9
```

***

1. Check that `fastqc` has been installed correctly by prompting the help message.
2. Check where `fastqc` was installed (hint: use `whereis`).

***

# YAML file {#yaml}
***
We have seen how simple it is to install tools using the `conda install` command.

In reality, we will want to install multiple packages at once for an analysis and create a clean environment for the packages. This can be simplified using a `.yml` file. The strucutre of a `.yml` file is:

1. Name: The name of the environment to be created.
2. Channels: Specifiy which channels conda should search when attempting to install packages.
3. Dependencies: Which packages you want to install (supports pinned version numbers).

***

In this weeks tutorial we want to create an environment for the quality control of sequencing reads. We will need `fastqc` and `multiqc` to generate HTML reports of sequencing statistics and a tool to perform adapter trimming and read filtering. Choosing a trimming tool is highly subjective, however I like the flexibility of `bbduk`, part of the `bbtools` suite.


Below I have provided an incomplete `.yml` file, where the name and channels have been specified for you. Please complete the `.yml` file by specifing the three tools referenced above.

```
name: QC
channels:
  - agbiome
  - bioconda
  - conda-forge
  - defaults
dependencies:
  -
  -
  -
```

***

Once you have filled out the `.yml` file, save it as `week1.yml`.

To create a conda environment using the `.yml` file, run the following command in the terminal:

```bash
conda env create -f week1.yml && conda clean -a
```

Conda should install the three tools under the environment `QC`.

### To Do:
> 1. Activate the environment.
> 2. Check all 3 tools have been installed correctly.
> 3. Print the path of the environments bin.
> 4. Export the environment using `conda env export > QC.yml`.
