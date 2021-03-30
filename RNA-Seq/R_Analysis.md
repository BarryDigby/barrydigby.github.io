---
title: R Analysis
layout: page
permalink: /RNA-Seq/R_Analysis
---

Before beginning the analysis in R, import the kallisto quantification directories to your local machine.

The `quant/` directory should have the following tree structure:
```bash
quant/
├── A375_1
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── A375_2
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── A375_3
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── A549_1
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── A549_2
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── A549_3
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── ctrl_1
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── ctrl_2
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
└── ctrl_3
    ├── abundance.h5
    ├── abundance.tsv
    └── run_info.json
```

`tarzip` the dirctory so we can move it to your `bactsrv` directory.

```
tar -zcvf quant.tar.gz quant/
```

Using `scp`, transfer the file to your own bactsrv account:

```
scp quant.tar.gz USERNAME@bactsrv.nuigalway.ie:/home/USERNAME
```

Close the connection to lugh and move to your local machine. Create a directory in your home directory and move there:

```
mkdir ~/RNA-Seq/
cd ~/RNA-Seq
```

Now download the tar file and unzip it.

```
scp USERNAME@bactsrv.nuigalway.ie:/home/USERNAME/quant.tar.gz ./ && tar -xvf quant.tar.gz
```

Create a sample table file. Save is as `samples.csv`

```
sample,condition,replicate
ctrl_1,control,1
ctrl_2,control,2
ctrl_3,control,3
A549_1,lung,1
A549_2,lung,2
A549_3,lung,3
A375_1,melanoma,1
A375_2,melanoma,2
A375_3,melanoma,3
```

***

# Rmd Document
Download the Rmd file here: [https://github.com/BarryDigby/barrydigby.github.io/blob/master/RNA-Seq/MA5512.Rmd](https://github.com/BarryDigby/barrydigby.github.io/blob/master/RNA-Seq/MA5112.Rmd).

or follow along with the Published HTML here: [https://rpubs.com/BarryDigby/747584](https://rpubs.com/BarryDigby/747584)
