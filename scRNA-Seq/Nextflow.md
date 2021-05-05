---
title: Nextflow
layout: page
permalink: /scRNA-Seq/Nextflow
---

# Walkthrough
The analysis using `kallisto` and `bustools` uses directory as outputs and inputs to downstream processes which we have not encountered before.

To call files in the directory, use parameter expansion: `${variable}`

Append a forward slash `/`, treating the variable as a directory: `${variable}/file.txt`

***

Save the below script as `test.nf` and run `nextflow run test.nf`

```bash
#!/usr/bin/env nextflow

/*
 * Use process to create output directory
 * with some files in it.
 */

process create_dir{

    output:
        file("bus_output") into dir_out

    script:
    """
    mkdir -p bus_output
    touch bus_output/output.bus
    touch bus_output/matrix.ec
    touch bus_output/transcripts.txt
    """
}

process stage_files{

    echo true

    input:
        file(bus_output) from dir_out

    output:
        stdout to out

    script:
    """
    printf "\n"
    printf "Ouput BUS file: ${bus_output}/output.bus \n"
    printf "Matrix file: ${bus_output}/matrix.ec \n"
    printf "Transcripts file: ${bus_output}/transcripts.txt \n"
    """
}
```
***

# Exercise
Complete the scRNA-Seq pipeline below. Refer to the `pipeline` section of this weeks tutorial to check process requirements.

This is a very staright forward pipeline as all variables are hard coded within the script.

Save the below script as `kallisto_bustools.nf` and run in your directory:

```bash
nextflow -bg \
run kallisto_bustools.nf \
-with-singularity /data/MSc/2020/MA5112/scRNA-Seq/container/scRNA.img
```

*N.B* I've tweaked the output directories slightly, all output files will end up in `bus_output/` which is slightly different from the interactive shell session we ran previously.

***

```bash
#!/usr/bin/env nextflow

reads = Channel.fromFilePairs("/data/MSc/2020/MA5112/scRNA-Seq/reads/*_R{1,2}_*")
index = Channel.value(file("/data/MSc/2020/MA5112/scRNA-Seq/reference/Homo_sapiens.cDNA.idx"))
whitelist = Channel.value(file("/data/MSc/2020/MA5112/scRNA-Seq/assets/10xv3_whitelist.txt"))
chemistry = "10xv3"
outdir = "."


process kallisto_bus{

    publishDir "${outdir}", mode:'copy'

    input:

    output:

    script:
    """
    kallisto bus \
    """
}

process bustools_correct{

    publishDir "${outdir}/bus_output", mode:'copy'

    input:

    output:

    script:
    """
    bustools correct \
    """
}


process bustools_sort{

    publishDir "${outdir}/bus_output", mode:'copy'

    input:

    output:

    script:
    """
    mkdir -p tmp

    bustools sort \
    """
}


process bustools_text{

    publishDir "${outdir}/bus_output", mode:'copy'

    input:

    output:

    script:
    """
    bustools text \
    """
}

```
