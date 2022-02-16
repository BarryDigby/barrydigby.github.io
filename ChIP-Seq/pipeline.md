# ChiP-Seq

### Housekeeping

Please do not run analyses in the `/home/` directory on lugh. All analyses must be run in the `/data/` directory.

Navigate to `/data/MSc/2022` and make a directory for yourself there using the naming convention `[first letter forename][surname]` eg: `bdigby`. 

If you have already created a directory like this with previous weeks tutorials, move it to this directory. Be sure to include the `-R` flag to add all sub-directories recursively. 

> Please delete all analysis files created in the `/home/` directory

## Tutorial files

All files required for the tutorial are present in `/data/MSc/2022/chip_scratch/files`:

* Reference Genome FASTA file
* Reference Genome GTF annotation file
* `SPT5_T0` samples (2 replicates)
* `SPT5_Input` controls (2 replicates)

> The `fastq` files hold ChiP-Seq data generated from _S.cerevisiae_ samples.  The files have been aggressively sub-sampled (100,000 reads per file) to speed up the analysis.

## Tutorial container

In addition to the analysis files, the container we are going to use is located in the container cache on lugh:

```bash
/data/containers/nfcore-chipseq-1.2.2.img
```

## Request resources

Before you start, request resources on a `MSC` compute node:

```bash
salloc -p MSC -n 1 -c 2
```

Check which node you have been assigned to using `squeue` and your username:

```bash
squeue -u bdigby
```

Shell into the compute node using `ssh`:

```bash
ssh compute[1/2/3]
```

**Be sure to change to your directory on `/data/MSc/2022/username` after this step.** 

Finally, shell into the container to access the suite of tools required for the tutorial:

```bash
module load singularity
singularity shell -B /data/ /data/containers/nfcore-chipseq-1.2.2.img
```

# Workflow

Copy the tutorial files to your directory on `/data/MSc/2022/username`:

```bash
cp -R /data/MSc/2022/chip_scratch/files/ .
```

# 1. BWA Index

Copy the `genome.fa` file to your current working directory and index the reference genome file using `bwa index`:

```bash
cp files/genome.fa .
bwa index -a bwtsw genome.fa
mkdir -p BWAIndex 
mv genome.fa* BWAIndex/
```


# 2. Quality Control

The `fastq` files provided for the analysis are already of high quality with no contamination. Regardless, we will generate `FastQC` reports for the sequencing reads as this is a step inherent in all genomic assays.

Make a directory for quality control files called `fastqc`, and run `fastqc` on the reads:

```bash
mkdir -p fastqc
fastqc files/*.fastq.gz --outdir fastqc/
```

# 3. Align Reads

We will align the sequencing reads to the reference genome using `bwa`. We require the sequencing files in pairs (both R1 and R2 are passed simultaneously), the reference genome FASTA file and the genome index files we previously generated. 



Recall there are 2 replicates for `SPT5_T0` and the input controls `SPT5_Input`. We will merge these replicates further downstream - it is not good practice to merge replicate fastq files prior to alignment (think about this critically). 



We will also need to specify the `@RG` header for `bwa` to specify membership in `bam` files i.e which replicates belong to `SPT5_T0` or `SPT5_Input` . This is crucial when merging `bam` files downstream.. 



### Symbolic Links

Create a directory called `bwa`. This is where we will perform alignment. 

We will use symbolic links to the reference index files and place them in the `bwa` directory:

```bash
mkdir -p bwa && cd bwa
ln -s ../BWAIndex/* .
```


Inspect the symbolic links, ask questions if what we are doing here is unclear. (I'm using symbolic links here as this is the underlying mechanism nextflow uses in its processes).

```bash
Singularity germline_vc.img:/data/MSc/2022/bdigby/chip_seq/bwa> ls -la
total 0
drwxr-xr-x  2 bdigby bdigby 155 Feb 15 15:01 .
drwxrwxr-x 10 bdigby bdigby 230 Feb 15 15:00 ..
lrwxrwxrwx  1 bdigby bdigby  21 Feb 15 15:01 genome.fa -> ../BWAIndex/genome.fa
lrwxrwxrwx  1 bdigby bdigby  25 Feb 15 15:01 genome.fa.amb -> ../BWAIndex/genome.fa.amb
lrwxrwxrwx  1 bdigby bdigby  25 Feb 15 15:01 genome.fa.ann -> ../BWAIndex/genome.fa.ann
lrwxrwxrwx  1 bdigby bdigby  25 Feb 15 15:01 genome.fa.bwt -> ../BWAIndex/genome.fa.bwt
lrwxrwxrwx  1 bdigby bdigby  25 Feb 15 15:01 genome.fa.pac -> ../BWAIndex/genome.fa.pac
lrwxrwxrwx  1 bdigby bdigby  24 Feb 15 15:01 genome.fa.sa -> ../BWAIndex/genome.fa.sa
```


### @RG Headers

> Please take your time with this, downstream processes will fail if you make an error here. 

Using `SPT5_T0_1_R1.fastq.gz` and `SPT5_T0_1_R2.fastq.gz` as an example, the corresponding `@RG` header will be:

```bash
@RG\tID:SPT5_T0_1\tSM:SPT5_T0\tPL:ILLUMINA\tLB:SPT5_T0_1\tPU:1
```

Notice the `SM:` sample identifier specifies this sample belongs to `SPT5_T0` of which there are two replicates. 

> The `ID` provided for each sample must be unique

## bwa mem

Align the files using the following command:

```bash
bwa mem -t 2 -M -R '@RG\tID:SPT5_T0_1\tSM:SPT5_T0\tPL:ILLUMINA\tLB:SPT5_T0_1\tPU:1' genome.fa ../files/SPT5_T0_1_R1.fastq.gz ../files/SPT5_T0_1_R2.fastq.gz | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o SPT5_T0_1.bam
```

> `genome.fa` is pointing to both the reference FASTA genome and the index files `genome.fa.*`

Repeat the above code for the 3 other samples `SPT5_T0_2_*`, `SPT5_Input_1_*` and `SPT5_Input_2_*`. Update the input `@RG` header, the input reads and output `bam` name for each iteration.

When you are finished the following `bam` files should be in your `bwa` directory:

```bash
SPT5_Input_1.bam  SPT5_Input_2.bam  SPT5_T0_1.bam  SPT5_T0_2.bam
```

## Automating the process

Automating this process is slightly tricky but solution is broken down below:

### Capturing read pairs

```bash
for file in ../files/*_R1.fastq.gz; do R1=$file; R2=${file/_R1.fastq.gz/}_R2.fastq.gz; echo $R1 $R2; done
```

gives:

```bash
../files/SPT5_Input_1_R1.fastq.gz ../files/SPT5_Input_1_R2.fastq.gz
../files/SPT5_Input_2_R1.fastq.gz ../files/SPT5_Input_2_R2.fastq.gz
../files/SPT5_T0_1_R1.fastq.gz ../files/SPT5_T0_1_R2.fastq.gz
../files/SPT5_T0_2_R1.fastq.gz ../files/SPT5_T0_2_R2.fastq.gz
```

which is exactly the input we want to pass to `bwa`, read pairs passed at the same time. 

The next step is to extract the ID name and Sample name for the `@RG` header. Build on your previous code (Indentation optional, helps readability..):

```bash
for file in ../files/*_R1.fastq.gz; do \
    
    R1=$file; \
    R2=${file/_R1.fastq.gz/}_R2.fastq.gz; \
    id=$(basename $file _R1.fastq.gz); \
    sample=${id%_*}; \
    
    echo $R1 $R2 $id $sample; \
    
done
```

gives:

```bash
../files/SPT5_Input_1_R1.fastq.gz ../files/SPT5_Input_1_R2.fastq.gz SPT5_Input_1 SPT5_Input
../files/SPT5_Input_2_R1.fastq.gz ../files/SPT5_Input_2_R2.fastq.gz SPT5_Input_2 SPT5_Input
../files/SPT5_T0_1_R1.fastq.gz ../files/SPT5_T0_1_R2.fastq.gz SPT5_T0_1 SPT5_T0
../files/SPT5_T0_2_R1.fastq.gz ../files/SPT5_T0_2_R2.fastq.gz SPT5_T0_2 SPT5_T0
```

We have everything we need for `bwa` automation, just wrap the for loop in the `bwa` command.

> Don't worry about the code used - the important part is that you are able to identify patterns in the file names and think critically about how they map to the `bwa` command. Once you know what you want to achieve, googling the code becomes much easier.

### Full Monty

```bash
for file in ../files/*_R1.fastq.gz; do \

    R1=$file; \
    R2=${file/_R1.fastq.gz/}_R2.fastq.gz; \
    id=$(basename $file _R1.fastq.gz); \
    sample=${id%_*}; \
    
    bwa mem -t 2 \
            -M \
            -R "@RG\tID:${id}\tSM:${sample}\tPL:ILLUMINA\tLB:${id}\tPU:1" \
            genome.fa \
            $R1 $R2 | samtools view -@ 2 -b -h -F 0x0100 -O BAM -o ${id}.bam; \

done
```

> You must place the `@RG` string in double quotes for the variable expansion to work.

#### Sanity Check

Make sure that the `@RG` headers are correctly formatted. You can use the command `samtools view -H <bam>` to inspect the header information of a `bam` file which contains chromosome information, `@RG` headers, and the history of commands used to generate the `bam` file. We will use `grep` to isolate the `@RG` headers:

```bash
for bam in *.bam; do samtools view -H $bam | grep "@RG"; done
```

> Note: If you recieve `bam` files for an analysis and want to know how they were produced, inspect the `@PG` lines in the header.

# 3. Samtools sort

After generating `bam` files for all of the samples the next step is to sort and index the files.

> run this in the `bwa/` directory containing the `bam` files generated in the previous step.

In the interest of time, use the for loop to do this step:

```bash
for bam in *.bam; do \

    prefix=$(basename $bam .bam);\
    
    samtools sort -@ 2 -o ${prefix}.sorted.bam $bam;\
    samtools index ${prefix}.sorted.bam
    
done
```

The `bam` files are now sorted by coordinates and indexed. This is done for computational efficiency and required for downstream steps.

# 4. picard

Now we can merge our replicates into one `bam` file before marking duplicates, the final preprocessing step prior to `macs2` analysis.

> At each step, we must index the output `bam` files

## 4a. MergeSamFiles

```bash
picard MergeSamFiles \
       INPUT=SPT5_T0_1.sorted.bam \
       INPUT=SPT5_T0_2.sorted.bam \
       OUTPUT=SPT5_T0.sorted.bam \
       SORT_ORDER=coordinate \
       VALIDATION_STRINGENCY=LENIENT
```

### Index 

```bash
samtools index SPT5_T0.sorted.bam
```

## 4b. MarkDuplicates

```bash
picard MarkDuplicates \
       INPUT=SPT5_T0.sorted.bam \
       OUTPUT=SPT5_T0.markdups.bam \
       ASSUME_SORTED=true \
       REMOVE_DUPLICATES=false \
       VALIDATION_STRINGENCY=LENIENT \
       METRICS_FILE=../fastqc/SPT5_T0.markdups.metrics.txt
```

### Index

```bash
samtools index SPT5_T0.markdups.bam
```

Repeat these 4 steps for the `SPT5_Input` samples. We will be using the `markdups.bam` files for peak calling.

This is seriously repetitive stuff..... . . . . . ..... . . .  . and we only have 4 replicates 0_o

# 5. Macs2

At long last.... 

move to the top of your directory tree for the analysis. make a new directory called `macs2` which is where we will store the results. 

```bash
mkdir -p macs2
```

run macs on our two bam files:

```bash
macs2 callpeak -t bwa/SPT5_T0.markdups.bam -c bwa/SPT5_Input.markdups.bam -f BAMPE --outdir macs2/ --name SPT5
```

# 6. Annotate Peaks

`annotatePeaks.pl` is a cool perl script that provides dense annotation of each peak using the input `GTF` file:

```bash
annotatePeaks.pl macs2/SPT5_summits.bed files/genome.fa -gid -gtf files/genes.gtf > SPT5_annotated_peaks.txt
```

Check your output file is the same:

```bash
wc -l SPT5_annotated_peaks.txt
714
```

```bash
head SPT5_annotated_peaks.txt
PeakID (cmd=annotatePeaks.pl SPT5_summits.bed ../files/genome.fa -gid -gtf ../files/genes.gtf)	Chr	Start	End	Strand	Peak Score	Focus Ratio/Region SizeAnnotation	Detailed Annotation	Distance to TSS	Nearest PromoterID	Entrez ID	Nearest Unigene	Nearest Refseq	Nearest Ensembl	Gene Name	Gene Alias	Gene Description	Gene Type
SPT5_peak_460	XII	460510	460510	+	838.196	NA	TTS (RDN25-2)	rRNA-TTS (RDN25-2)	834	RDN5-1	RDN5-1	RDN5-1			RDN5-1		rRNA
SPT5_peak_465	XII	490122	490122	+	218.423	NA	promoter-TSS (YLR162W-A)	protein_coding-promoter-TSS (YLR162W-A)	-284	YLR162W-A	YLR162W-A	YLR162W-A			RRT15			protein_coding
SPT5_peak_510	XII	1071609	1071609	+	139.342	NA	promoter-TSS (YLR466C-A)	protein_coding-promoter-TSS (YLR466C-A)	101	YLR466C-B	YLR466C-B	YLR466C-B			YLR466C-B			protein_coding
SPT5_peak_662	XV	1091237	1091237	+	134.376	NA	promoter-TSS (YOR396C-A)	protein_coding-promoter-TSS (YOR396C-A)	-740	YOR396C-A	YOR396C-A	YOR396C-A			YOR396C-A			protein_coding
SPT5_peak_332	VIII	68	68	+	80.9434	NA	promoter-TSS (YHL050W-A)	protein_coding-promoter-TSS (YHL050W-A)	-743	YHL050W-A	YHL050W-A	YHL050W-A			YHL050W-A			protein_coding
SPT5_peak_445	XII	5739	5739	+	77.7621	NA	promoter-TSS (YLL066W-A)	protein_coding-promoter-TSS (YLL066W-A)	134	YLL066W-B	YLL066W-B	YLL066W-B			YLL066W-B			protein_coding
SPT5_peak_664	XVI	6457	6457	+	51.9798	NA	promoter-TSS (YPL283C)	protein_coding-promoter-TSS (YPL283C)	-450	YPL283C	YPL283C	YPL283C		YRF1-7			protein_coding
SPT5_peak_446	XII	11181	11181	+	51.3152	NA	promoter-TSS (YLL065W)	protein_coding-promoter-TSS (YLL065W)	-545	YLL065W	YLL065W	YLL065W		YLL065W			protein_coding
SPT5_peak_331	VII	1084549	1084549	+	48.2078	NA	promoter-TSS (YGR296W)	protein_coding-promoter-TSS (YGR296W)	-315	YGR296W	YGR296W	YGR296W		YRF1-3			protein_coding
```