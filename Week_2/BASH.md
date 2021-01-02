---
title: Exome Sequencing Analysis
layout: page
permalink: /Week_2/BASH
---

A detailed workflow has been provided for you to run a germline variant calling analysis on LUGH before attempting to create a nextflow script.

1. [Genome Indexing](#index)
    - [BWA Index](#bwaidx)
    - [Samtools Index](#faidx)
    - [Picard CreateSequenceDictionary](#seqdict)
2. [Align Reads](#align)
3. [Mark Duplicates](#markdup)
4. [BQSR](#bqsr)
    - [Base Recalibration](#baserecal)
    - [Apply BQSR](#applybqsr)
5. [Germline Variant Calling](#germline_vc)
    - [Haplotype Caller](#haplotype)
    - [Genotype VCFs](#genotype)
6. [Subset Variants](#subset)
    - [SNPs](#subsetsnp)
    - [INDELs](#subsetindel)
7. [Filter Variants](#filter)
    - [SNPs](#filtersnp)
    - [INDELs](#filterindel)
8. [Merge VCFs](#mergevcf)
9. [Annotate Variants](#annotate)

***

## Compute Resources
**Do not run this analysis on the head node**.

Request computing resources on LUGH:

```bash
salloc -p MSC -c 1 -n 1
```

```bash
ssh compute0{1..3}
```

***

Tools required for the analysis are available in a pre-prepared container for the tutorial. Invoke an interactive shell session within the container.

```
singularity shell -B /data /path/to/container
```

Ask for help if you are uncertain if you have completed this step correctly.

# 1. Genome Index {#index}
**Due to time constraints all indexing has been performed for you. Skip to step 2.**  

***

## **BWA Index** {#bwaidx}
BWA requires building an index for your reference genome to allow computationally efficient searches of the genome during sequence alignment.

##### *Inputs*
- Reference Genome

```bash
bwa index -a bwtsw GRCh38.fasta
```

##### *Flags*
* `-a`: Construction algorithm (bwtsw, is or rb2). For large genomes, specify `-a bwtsw`. If you are unsure, omit this flag and `bwa` will determine the correct algorithm to use.

##### *Outputs*
BWA indexing produces 5 files with the file extensions `*.amb`, `*.ann`, `*btw`, `*.pac` & `*.sa`.

***

## **Samtools Index** {#faidx}
Indexes (or queries) the reference sequence.

##### *Inputs*
- Reference Genome

```bash
samtools faidx GRCh37.fasta
```
##### *Outputs*
Samtools fasta indexing generates a `*.fai` file.

***

## **Picard CreateSequenceDictionary** {#seqdict}
Creates a sequence dictionary of the reference fasta file specifying the chromosomes and chromosome sizes. Required for downstream analysis tools by GATK.

##### *Inputs*
- Reference Genome

```bash
picard CreateSequenceDictionary \
       R=GRCh37.fasta \
       O=GRCh37.dict
```

##### *Outputs*
A dictionary file. You may name this whatever you want, however common convention dictates `*.dict`.

# 2. **Align Reads** {#align}
Map the reads to the reference genome using `bwa mem`. The SAM files produced have reads in the order that the sequences occurred in the input FASTQ files i.e in randomn order. `samtools sort` orders the reads by their leftmost coordinate i.e in 'genome order'. This is a requirement for downstream tools.

We will also convert the SAM file to to its binary counterpart: a BAM file.

##### *Inputs*
* Reference Genome
* FASTQ Reads

```bash
bwa mem \
    -K 10000000 \
    -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' \
    -t 8 \
    GRCh37.fasta \
    ../reads/subsample_r1.fastq.gz \
    ../reads/subsample_r2.fastq.gz \
    | samtools sort - > subsample.bam
```

##### *Flags*
* `-K`: process INT input bases in each batch regardless of nThreads (for reproducibility)
* `-R`: Read Group identifier. This information is key for downstream GATK functionality, GATK will not work without a read group tag.
* `-t`: Number of threads

##### *Outputs*
The script passes the output of `bwa mem` (`subsample.sam`) directly to `samtools sort` by using the pipe command (`|`). A dash `-` is used to indicate the incoming file from the previous process for samtools. When working with full sized samples, allocate more threads to samtools like so `samtools sort --threads 8 - > subsample.bam`.

*N.B*: `bwa mem` requires all `bwa idx` output files & the reference genome to be present in the working directory for the tool to run.

***

# 3. Mark Duplicates {#markdup}
*"Almost all statistical models for variant calling assume some sort of independence between measurements. The duplicates (if one assumes that they arise from PCR artifact) are not independent. This lack of independence will usually lead to a breakdown of the statistical model and measures of statistical significance that are incorrect"* -- Sean Davis.

##### *Inputs*
- Sorted BAM file.

```bash
gatk --java-options -Xmx2g \
     MarkDuplicates \
     --MAX_RECORDS_IN_RAM 50000 \
     --INPUT subsample.bam \
     --METRICS_FILE subsample.bam.metrics \
     --TMP_DIR ./ \
     --ASSUME_SORT_ORDER coordinate \
     --CREATE_INDEX true \
     --OUTPUT subsample.markdup.bam
```

##### *Outputs*
- `--METRICS_FILE`: Statistics of duplication events in sequencing reads.
- `--OUTPUT`: A BAM file with collapsed duplicates.

# 4. BQSR {#bqsr}
Base Quality Score Recalibration. A data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call.

## BaseRecalibrator {#baserecal}
To build the recalibration model, `BaseRecalibrator` goes through all of the reads in the input BAM file and tabulates data about the following features of the bases:

- read group the read belongs to
- quality score reported by the machine
- machine cycle producing this base (Nth cycle = Nth base from the start of the read)
- current base + previous base (dinucleotide)

For each bin, the tool counts the number of bases within the bin and how often such bases mismatch the reference base, excluding loci known to vary in the population, according to the known variants resource (typically dbSNP, Mills_KG_gold). This information is output to a recalibration file in GATKReport format.

***

##### *Inputs*
- Marked Duplicates BAM file
- Exome Interval list
- Reference Genome
- Known variants (SNPs + INDELs)

```bash
gatk --java-options -Xmx2g \
     BaseRecalibrator \
     -I subsample.markdup.bam \
     -O subsample.recal.table \
     -L ../assets/exome.bed.interval_list \
     --tmp-dir . \
     -R GRCh37.fa \
     --known-sites ../assets/Mills_KG_gold.indels.b37.vcf.gz \
     --known-sites ../assets/dbsnp_138.b37.vcf.gz
```

##### *Outputs*
- `-O`: Recalibration table file. Required for next step.

***

## ApplyBQSR {#applybqsr}
`ApplyBQSR` goes through all the reads again, using the recalibration table file to adjust each base's score based on which bins it falls in. So effectively the new quality score is:

- The sum of the global difference between reported quality scores and the empirical quality
- Plus the quality bin specific shift
- Plus the cycle x qual and dinucleotide x qual effect

Following recalibration, the read quality scores are much closer to their empirical scores than before. This means they can be used in a statistically robust manner for downstream processing, such as variant calling. In addition, by accounting for quality changes by cycle and sequence context, we can identify truly high quality bases in the reads, often finding a subset of bases that are Q30 even when no bases were originally labeled as such.

##### *Input*
- Marked Duplicates BAM file
- Exome Interval list
- Reference Genome
- Recalibration table file (previous step).

```bash
gatk --java-options -Xmx2g \
     ApplyBQSR \
     -I subsample.markdup.bam \
     -O subsample.recal.bam \
     -L ../assets/exome.bed.interval_list \
     -R GRCh37.fasta \
     --bqsr-recal-file subsample.recal.table
```

##### *Outputs*
- `-O`: Recalibrated, Marked Duplicate BAM file.

***

Tidy up the indexed BAM file & retrieve the statistics of the recalibrated BAM file.

```bash
mv subsample.recal.bai subsample.recal.bam.bai
```

```bash
samtools stats subsample.recal.bam > subsample.recal.stats
```

# 5. Germline Variant Calling {#germline_vc}
Identify germline short variants (SNPs and Indels) in an individual, or in a cohort. This tutorial is focused on a single sample germline variant calling analysis. Not to be confused with somatic variant calling which uses matched tumour - normal samples, requiring a different workflow.

Performing a single sample analysis and a joint cohort analysis follow the same steps covered below.

## Haplotype Caller {#haplotype}
Call germline SNPs and indels via local re-assembly of haplotypes via:

1. *Defining active regions*
    - The program determines which regions of the genome it needs to operate on (active regions), based on the presence of evidence for variation.
2. *Determine haplotypes by assembly of the active region*
    - For each active region, the program builds a De Bruijn-like graph to reassemble the active region and identifies what are the possible haplotypes present in the data. The program then realigns each haplotype against the reference haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites.
3. *Determine likelihoods of the haplotypes given the read data*
    - For each active region, the program performs a pairwise alignment of each read against each haplotype using the PairHMM algorithm. This produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.
4. *Assign sample genotypes*
    - For each potential variant site, the program applies Bayes' rule, using the likelihoods of alleles given the read data to calculate the likelihoods of each genotype per sample given the read data observed for that sample. The most likely genotype is then assigned to the sample.

***

##### *Inputs*
- Reference Genome
- Recalibrated BAM
- Exome Interval list
- Known SNPs (dbSNP)

```bash
gatk --java-options -Xmx2g \
     HaplotypeCaller \
     -R GRCh37.fasta \
     -I subsample.recal.bam \
     -L ../assets/exome.bed.interval_list \
     -D ../assets/dbsnp_138.b37.vcf.gz \
     -O subsample.g.vcf \
     -ERC GVCF
```

##### *Flags*
- `ERC`: Mode for emitting reference confidence scores. The key difference between a regular VCF and a GVCF is that the GVCF has records for all sites, whether there is a variant call there or not. The goal is to have every site represented in the file in order to do joint analysis of a cohort in subsequent steps.

##### *Outputs*
- `-O`: GVCF as described above.

***

## GenotypeVCFs {#genotype}
`GenotypeVCFs` is designed to perform joint genotyping on a single input, which may contain one or many samples. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.

##### *Inputs*
- Reference Genome
- Exome Interval list
- GVCF
- Known SNPs (dbSNP)

```bash
gatk --java-options -Xmx2g \
     GenotypeGVCFs \
     -R GRCh37.fasta \
     -L ../assets/exome.bed.interval_list \
     -D ../assets/dbsnp_138.b37.vcf.gz \
     -V subsample.g.vcf \
     -O subsample.vcf
```

##### *Outputs*
- `-O`: VCF file containing variants.

***

# 6. Subset Variants {#subset}
Before applying filtering thresholds to the called variants, subset the VCF file for both SNPs and INDELs using `SelectVariants`.

## SNPs {#subsetsnp}

##### *Inputs*
- Reference Genome
- VCF file produced by `GenotypeVCFs`

```bash
gatk SelectVariants \
     -R GRCh37.fasta \
     -V subsample.vcf \
     -O subsample.snps.vcf.gz \
     -select-type SNP
```

##### *Outputs*
- `-O`: VCF file containing SNPs

## INDELs {#subsetindel}

##### *Inputs*
- Reference Genome
- VCF file produced by `GenotypeVCFs`

```bash
gatk SelectVariants \
     -R GRCh37.fasta \
     -V subsample.vcf \
     -O subsample.indels.vcf.gz \
     -select-type INDEL
```

##### *Outputs*
- `-O`: VCF file containing INDELs

# 7. Filter Variants {#filter}
Apply filtering to selected variants generated in the previous step.

## SNPs {#filtersnp}
Please refer to the following [blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) by GATK regarding hard-filtering SNPs.

##### *Inputs*
- Reference Genome
- VCF file containing SNPs

```bash
gatk VariantFiltration \
     -R GRCh37.fasta \
     -V subsample.snps.vcf.gz \
     -O subsample.filt_snps.vcf \
     --filter-expression "QD < 2.0" \
     --filter-name "filterQD_lt2.0" \
     --filter-expression "MQ < 25.0" \
     --filter-name "filterMQ_lt25.0" \
     --filter-expression "SOR > 3.0" \
     --filter-name "filterSOR_gt3.0" \
     --filter-expression "MQRankSum < -12.5" \
     --filter-name "filterMQRankSum_lt-12.5" \
     --filter-expression "ReadPosRankSum < -8.0" \
     --filter-name "filterReadPosRankSum_lt-8.0"
```

##### *Flags*
- `--filter-expression`: Filtering to be applied to SNPs.
- `--filter-name`: Annotation of SNPs failing filtering threshold written to VCF file.

##### *Outputs*
- `-O`: VCF file with annotated SNPs failing the selected filtering thresholds. SNPs that pass all thrtesholds will be marked with `PASS`.

***

## INDELs {#filterindel}
Please refer to the following [blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) by GATK regarding hard-filtering INDELs.

##### *Inputs*
- Reference Genome
- VCF file containing INDELs

```bash
gatk VariantFiltration \
     -R GRCh37.fasta \
     -V subsample.vcf \
     -O subsample.filt_indels.vcf \
     --filter-expression "QD < 2.0" \
     --filter-name "filterQD" \
     --filter-expression "SOR > 10.0" \
     --filter-name "filterSOR_gt10.0" \
     --filter-expression "ReadPosRankSum < -20.0" \
     --filter-name "filterReadPosRankSum"
```
##### *Flags*
- `--filter-expression`: Filtering to be applied to INDELs.
- `--filter-name`: Annotation of INDELs failing filtering threshold written to VCF file.

##### *Outputs*
- `-O`: VCF file with annotated INDELs failing the selected filtering thresholds. INDELs that pass all thrtesholds will be marked with `PASS`.

# 8. Merge VCFs {#mergevcf}
`MergeVCFs` combines multiple variant files into a single variant file.

##### *Inputs*
- Filtered SNPs VCF file
- Filtered INDELs VCF file

```bash
gatk MergeVcfs \
     -I subsample.filt_indels.vcf \
     -I subsample.filt_snps.vcf \
     -O filtered_sample.vcf
```

##### *Outputs*
- `-O`: The same VCF file output by `GenotypeVCFs`, with added annotations denoting `PASS` or specifying filtering thresholds failed.

# 9. Annotate Variants {#annotate}
Run in your own time, this step takes a long time. Usually we submit this job to SLURM with much more resources than we have initially requested.

```bash
snpEff GRCh37.75 \
       -csvStats subsample.snpEff.csv \
       -nodownload \
       -dataDir /data/snpEff/ \
       -canon \
       -v filtered_sample.vcf > subsample.snpEff.ann.vcf
```
