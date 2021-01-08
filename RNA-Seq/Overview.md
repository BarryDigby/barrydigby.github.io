---
title: Overview
layout: page
permalink: /RNA-Seq/Overview
---

Running a RNA-Seq quantification analysis on LUGH is extremely straight forward thanks to the simplicity of `kallisto`. The majority of this tutorial will be spent on downstream differential expression analysis in R.

1. Nextflow tutorial + exercise *~1 hour*
2. DESeq2 Differential Expression Analysis in R *~ 1 hour*

Please install the following packages in R Studio:

```
dplyr
biomaRt
tximport
rhdf5
gplots
DESeq2
apeglm
RColorBrewer
IHW
PCAtools
EnhancedVolcano
ComplexHeatmap
circlize
fgsea
tidyverse
ggpubr
```

Unfortunately we cannot mount a conda environment in R Studio (to my knowledge) thus each package must be installed manually.

***

# Dataset

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/rna-seq/rnaseq-data.png" width="100%" height="100%"/>
</center>

The researcher has captured exosomes from human malignant melanoma cell line (A375) and the human lung cancer cell line (A549) and added them to primary human epidermal melanocytes, (NHEM-c cells). There are triplicates of each cell line:

- NHEM-c (Control)
- NHEM-c + A549 exosomes (Lung)
- NHEM-c + A375 exosomes (Melanoma)

The researcher is interested in characterising the transcriptional profile of each transfected cell line, in an effort to elucidate the tumor pathways mediated by A375 & A549 tumor exosomes.
