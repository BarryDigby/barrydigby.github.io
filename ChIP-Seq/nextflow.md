# Prerequisites

```bash
module load singularity/3.4.1
module load java/1.8.0
```

The `nextflow` executable has been placed in `/data/MSc/2022/` for you. Add this to your bin:

```bash
cp /data/MSc/2022/nextflow ~/bin/
```

Double check you can access the executable:

```bash
nextflow -v
nextflow version 21.04.0-edge.5543
```

## Config file

We need to tell nextflow to run on the `MSC` queue. 

You should already have a config file present in your `~/.nextflow` directory (we used this to run the variant calling pipeline).

Replace the following line:

```bash
queue = { task.cpus > 8 ? 'highmem' : 'normal' }
```

with 

```bash
queue = 'MSC'
```

# Staging Files

The first step we will cover is reading files into `nextflow` scripts and the data structure they are stored in.

Copy the contents of the code below, save it in a script called `stage_files.nf` and run it:

`nextflow run stage_files.nf`

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

ch_reads.view()
ch_fasta.view()
ch_gtf.view()
```

# Building your first process

We will build on our previous script and include a process to perform indexing using `bwa`.

We will also add a `parameter` called `outdir` which we will use to specify the path to write files to. 

> parameters can be hardcoded in the script, passed via the command line, or config file.

The flag `-with-singularity` needs to be included in the command line, the script needs to know where to find the tools to run the process.

`nextflow run stage_files.nf --outdir "./" -with-singularity "/data/containers/nfcore-chipseq-1.2.2.img"`

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

process BWA_Index{
    tag "Indexing $fasta"
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex/*") into ch_index

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir -p BWAIndex
    mv ${fasta}* BWAIndex/
    """
}

ch_index.view()
```

> Note: generally you cannot reuse channels twice. Running `.view()` 'consumes' the channel, so after your sanity check remove the `.view()` line.

# Aligning Reads

Build on the previous script by adding alignment using `bwa`.

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

process BWA_Index{
    tag "Indexing $fasta"
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex/*") into ch_index

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir -p BWAIndex
    mv ${fasta}* BWAIndex/
    """
}


process BWA_Align{
    tag "Aligning $base to genome"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(reads) from ch_reads
    file(index) from ch_index.collect()

    output:
    tuple val(base), file("*.bam") into ch_aligned_bams

    script:
    RG = "\'@RG\\tID:${base}\\tSM:${base.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${base}\\tPU:1\'"
    genome = index[0]
    """
    bwa mem \
        -t 2 \
        -M \
        -R $RG \
        $genome \
        $reads | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o ${base}.bam
    """
}
```

# Samtools Sort 


Update the script to include the sorting of aligned `bam` files:

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

process BWA_Index{
    tag "Indexing $fasta"
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex/*") into ch_index

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir -p BWAIndex
    mv ${fasta}* BWAIndex/
    """
}


process BWA_Align{
    tag "Aligning $base to genome"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(reads) from ch_reads
    file(index) from ch_index.collect()

    output:
    tuple val(base), file("*.bam") into ch_aligned_bams

    script:
    RG = "\'@RG\\tID:${base}\\tSM:${base.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${base}\\tPU:1\'"
    genome = index[0]
    """
    bwa mem \
        -t 2 \
        -M \
        -R $RG \
        $genome \
        $reads | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o ${base}.bam
    """
}


process Samtools_Sort{
    tag "Sorting sample $base"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_aligned_bams

    output:
    tuple val(base), file("*.{bam,bai}") into ch_aligned_bams_sorted

    script:
    """
    samtools sort -@ 2 -o ${base}.sorted.bam $bam
    samtools index ${base}.sorted.bam
    """
}
```

# Merging BAMs

One of the hardest things to conceptualise in nextflow is merging tuples. The schematic below shows how to merge our sorted `bam` files according to a common grouping key i.e the `Sample` identifier, which we can grab from the file names. 

There exists multiple ways to go about this!

> If you need help with this sort of operation in your final assignment please do not hesitate to ask me for help.

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/mergesams.png" width="100%" height="100%"/>
</center>

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

process BWA_Index{
    tag "Indexing $fasta"
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex/*") into ch_index

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir -p BWAIndex
    mv ${fasta}* BWAIndex/
    """
}


process BWA_Align{
    tag "Aligning $base to genome"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(reads) from ch_reads
    file(index) from ch_index.collect()

    output:
    tuple val(base), file("*.bam") into ch_aligned_bams

    script:
    RG = "\'@RG\\tID:${base}\\tSM:${base.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${base}\\tPU:1\'"
    genome = index[0]
    """
    bwa mem \
        -t 2 \
        -M \
        -R $RG \
        $genome \
        $reads | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o ${base}.bam
    """
}


process Samtools_Sort{
    tag "Sorting sample $base"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_aligned_bams

    output:
    tuple val(base), file("*.{bam,bai}") into ch_aligned_bams_sorted

    script:
    """
    samtools sort -@ 2 -o ${base}.sorted.bam $bam
    samtools index ${base}.sorted.bam
    """
}


ch_aligned_bams_sorted
    .map{ it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
    .groupTuple(by: [0])            
    .map{ it -> [ it[0], it[1].flatten() ] }
    .set{ ch_merge_bams }


process Merge_Bams{
    tag "Merging $base BAM files"
    publishDir "${params.outdir}/merged_bams", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_merge_bams

    output:
    tuple val(base), file("*.{bam,bai}") into ch_merged_bams

    script:
    bams = bam.findAll{ it.toString().endsWith('.bam') }.sort()
    """
    picard MergeSamFiles \
           ${'INPUT='+bams.join(' INPUT=')} \
           OUTPUT=${base}.sorted.bam \
           SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT \
           
    samtools index ${base}.sorted.bam
    """
}
```

# MarkDuplicates

This process is straight forward now that we have merged our `bam` files. 

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

process BWA_Index{
    tag "Indexing $fasta"
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex/*") into ch_index

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir -p BWAIndex
    mv ${fasta}* BWAIndex/
    """
}


process BWA_Align{
    tag "Aligning $base to genome"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(reads) from ch_reads
    file(index) from ch_index.collect()

    output:
    tuple val(base), file("*.bam") into ch_aligned_bams

    script:
    RG = "\'@RG\\tID:${base}\\tSM:${base.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${base}\\tPU:1\'"
    genome = index[0]
    """
    bwa mem \
        -t 2 \
        -M \
        -R $RG \
        $genome \
        $reads | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o ${base}.bam
    """
}


process Samtools_Sort{
    tag "Sorting sample $base"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_aligned_bams

    output:
    tuple val(base), file("*.{bam,bai}") into ch_aligned_bams_sorted

    script:
    """
    samtools sort -@ 2 -o ${base}.sorted.bam $bam
    samtools index ${base}.sorted.bam
    """
}


ch_aligned_bams_sorted
    .map{ it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
    .groupTuple(by: [0])            
    .map{ it -> [ it[0], it[1].flatten() ] }
    .set{ ch_merge_bams }


process Merge_Bams{
    tag "Merging $base BAM files"
    publishDir "${params.outdir}/merged_bams", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_merge_bams

    output:
    tuple val(base), file("*.bam"), file("*.bai") into ch_merged_bams

    script:
    bams = bam.findAll{ it.toString().endsWith('.bam') }.sort()
    """
    picard MergeSamFiles \
           ${'INPUT='+bams.join(' INPUT=')} \
           OUTPUT=${base}.sorted.bam \
           SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT \
           
    samtools index ${base}.sorted.bam
    """
}


process MarkDuplicates{
    tag "Marking Duplicates in $base"
    publishDir "${params.outdir}/markdups", mode:'copy'

    input:
    tuple val(base), file(bam), file(bai) from ch_merged_bams

    output:
    tuple val(base), file("*.bam"), file("*.bai") into ch_markdup_bams

    script:
    """
    picard MarkDuplicates \
           INPUT=$bam \
           OUTPUT=${base}.markdups.bam \
           ASSUME_SORTED=true \
           REMOVE_DUPLICATES=false \
           VALIDATION_STRINGENCY=LENIENT \
           METRICS_FILE=${base}.markdups.metrics.txt

    samtools index ${base}.markdups.bam
    """
}
```

# Macs2

I was unable to expicitly pass `SPT5_T0` and `SPT5_Input` to the `-t` and `-c` flags dynamically.

This will work for the tutorial but is not scalable to other samples.. 

```groovy
#!/usr/bin/env nextflow

ch_reads = Channel.fromFilePairs("/data/MSc/2022/chip_scratch/files/*_R{1,2}.fastq.gz")
ch_fasta = Channel.value(file("/data/MSc/2022/chip_scratch/files/genome.fa"))
ch_gtf = Channel.value(file("/data/MSc/2022/chip_scratch/files/genes.gtf"))

process BWA_Index{
    tag "Indexing $fasta"
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex/*") into ch_index

    script:
    """
    bwa index -a bwtsw $fasta
    mkdir -p BWAIndex
    mv ${fasta}* BWAIndex/
    """
}


process BWA_Align{
    tag "Aligning $base to genome"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(reads) from ch_reads
    file(index) from ch_index.collect()

    output:
    tuple val(base), file("*.bam") into ch_aligned_bams

    script:
    RG = "\'@RG\\tID:${base}\\tSM:${base.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${base}\\tPU:1\'"
    genome = index[0]
    """
    bwa mem \
        -t 2 \
        -M \
        -R $RG \
        $genome \
        $reads | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o ${base}.bam
    """
}


process Samtools_Sort{
    tag "Sorting sample $base"
    publishDir "${params.outdir}/bwa", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_aligned_bams

    output:
    tuple val(base), file("*.{bam,bai}") into ch_aligned_bams_sorted

    script:
    """
    samtools sort -@ 2 -o ${base}.sorted.bam $bam
    samtools index ${base}.sorted.bam
    """
}


ch_aligned_bams_sorted
    .map{ it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
    .groupTuple(by: [0])            
    .map{ it -> [ it[0], it[1].flatten() ] }
    .set{ ch_merge_bams }


process Merge_Bams{
    tag "Merging $base BAM files"
    publishDir "${params.outdir}/merged_bams", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_merge_bams

    output:
    tuple val(base), file("*.bam"), file("*.bai") into ch_merged_bams

    script:
    bams = bam.findAll{ it.toString().endsWith('.bam') }.sort()
    """
    picard MergeSamFiles \
           ${'INPUT='+bams.join(' INPUT=')} \
           OUTPUT=${base}.sorted.bam \
           SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT \
           
    samtools index ${base}.sorted.bam
    """
}


process MarkDuplicates{
    tag "Marking Duplicates in $base"
    publishDir "${params.outdir}/markdups", mode:'copy'

    input:
    tuple val(base), file(bam), file(bai) from ch_merged_bams

    output:
    tuple val(base), file("*.bam"), file("*.bai") into ch_markdup_bams

    script:
    """
    picard MarkDuplicates \
           INPUT=$bam \
           OUTPUT=${base}.markdups.bam \
           ASSUME_SORTED=true \
           REMOVE_DUPLICATES=false \
           VALIDATION_STRINGENCY=LENIENT \
           METRICS_FILE=${base}.markdups.metrics.txt

    samtools index ${base}.markdups.bam
    """
}


ch_markdup_bams
    .map{ it -> [ it[0].split('_')[0], it[1] ] }
    .groupTuple( by:[0] )
    .map{ it -> [ it[0], it[1].flatten() ] }
    .set{ ch_macs_bams }


process Macs2{
    tag "Performing Peak Calling on $base"
    publishDir "${params.outdir}/macs2", mode:'copy'

    input:
    tuple val(base), file(bam) from ch_macs_bams

    output:
    file("${base}_*") into ch_macs2_output
    file("${base}_summits.bed") into ch_macs2_bed

    script:
    """
    macs2 callpeak \
          -t SPT5_T0.markdups.bam \
          -c SPT5_Input.markdups.bam \
          -f BAMPE \
          --outdir . \
          --name $base
    """
}
```

# Assignment

Finish the script by adding the final process to annotate the peaks using `annotatePeaks.pl`. 

I have placed the `SPT5_summits.bed` file in it's own channel - calling it should be straightforward. 
