---
title: R Analysis
layout: page
permalink: /scRNA-Seq/R_Analysis
---

Before beginning the analysis in R, we have to convert the kallisto output files into a format compatible with Seurat.

Your results directory `analysis/` should have the following tree structure:

```bash
analysis/
├── counts
│   └── output.sort.txt
├── raw
│   ├── bus_output
│   │   ├── matrix.ec
│   │   ├── output.bus
│   │   ├── run_info.json
│   │   └── transcripts.txt
│   └── kallisto.log
└── sorted
    └── bus_output
        ├── matrix.ec
        ├── output.bus
        ├── output.corrected.bus
        ├── output.sort.bus
        ├── run_info.json
        └── transcripts.txt
```

We are going to use a python script `make_mtx.py`, available at the following  [repository](https://github.com/BarryDigby/barrydigby.github.io/tree/master/scRNA-Seq)

Run `python make_mtx.py` to generate a seurat directory containing the files needed for the analysis.

`tarzip` the dirctory so we can move it to your `bactsrv` directory.

```
tar -zcvf seurat.tar.gz seurat/
```

Using `scp`, transfer the file to your own bactsrv account:

```
scp seurat.tar.gz USERNAME@bactsrv.nuigalway.ie:/home/USERNAME
```

Close the connection to lugh and move to your local machine. Create a directory in your home directory and move there:

```
mkdir ~/scRNA-Seq/
cd ~/scRNA-Seq
```

Now download the tar file and unzip it.

```
scp USERNAME@bactsrv.nuigalway.ie:/home/USERNAME/seurat.tar.gz ./ && tar -xvf seurat.tar.gz
```

***

The R markdown document for the analysis is available at the previous [repository](https://github.com/BarryDigby/barrydigby.github.io/tree/master/scRNA-Seq) link. 
