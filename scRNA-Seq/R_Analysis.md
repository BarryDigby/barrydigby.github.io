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
    ├── output.corrected.bus
    └── output.sort.bus
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

*N.B* one of the ENSEMBL gene identifiers has no gene symbol, and will create an error in seurat.

Identify the row and remove it:
```
awk 'NF!=2' genes.tsv
ENSG00000237235.2
grep -v "ENSG00000237235.2" genes.tsv > tmp.tsv && rm genes.tsv && mv tmp.tsv genes.tsv
```

Look up the gene ID (`ENSG00000237235.2`) in google and add back to the `genes.tsv` file using a text editor, using a tab delimiter.

Realistically we would do this in biomaRt to automate the process...
