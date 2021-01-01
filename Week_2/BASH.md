---
title:
layout: page
permalink: /Week_2/BASH
---

# Exome Sequencing Analysis
We are going to analyse a subset of a patients whole exome sequencing data to identify variants in protein coding genes, attempting to elucidate possible genetic diseases.

# Workflow
Request resources on the MSC queue. **Do not run this analysis on the head node** - unless you want an email from Declan :] !

```bash
salloc -p MSC -c 1 -n 1
```

Next, start an interactive shell with the container provided for the analysis.

```
singularity shell -B /data /path/to/container
```

### 1. Genome Index
As with any analysis, the reference genome must be indexed to quickly extract alignments overlapping particular genomic regions.

**Due to time constraints indexing has been performed for you. Skip to step 2.**

### BWA Index
```bash
bwa index -a bwtsw GRCh38.fasta
```

* `-a`: Construction algorithm (bwtsw, is or rb2). For large genomes, specify `-a bwtsw`. If you are unsure, omit this flag and `bwa` will determine the correct algorithm to use.

##### Output
Indexing using `bwa` produces 5 files: `*.amb`, `*.ann`, `*btw`, `*.pac` & `*.sa`.

***

### Samtools Index
```bash
samtools faidx GRCh37.fasta
```
##### Output
Samtools generates a `*.fai` file (<strong>fa</strong>sta <strong>i</strong>ndexed).

***

### Picard CreateSequenceDictionary
```bash
picard CreateSequenceDictionary R=GRCh37.fasta O=GRCh37.dict
```

***

### 2. Map Reads
Map the reads to the reference genome using `bwa mem`, and pipe the output to `samtools sort`.

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

* `-K`: process INT input bases in each batch regardless of nThreads (for reproducibility)
* `-R`: Read Group identifier. This information is key for downstream GATK functionality, GATK will not work without a read group tag.
* `-t`: Number of threads

The script above passes the output of `bwa mem` (`subsample.sam`) directly to `samtools sort` by using the pipe command (`|`). Samtools uses a dash `-` to indicate the incoming file from the previous process. When working with full sized samples, allocate more threads to samtools like so `samtools sort --threads 8 - > subsample.bam`.

**NOTE** `bwa mem` requires all genome index (both `bwa index`, `samtools index`) files to be present in the working directory for the tool to run. This will be important when designing the nextflow script.


### 3. Mark Duplicates

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

### 4. BQSR

### BaseRecalibrator
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


### ApplyBQSR
```bash
gatk --java-options -Xmx2g \
     ApplyBQSR \
     -I subsample.markdup.bam \
     -O subsample.recal.bam \
     -L ../assets/exome.bed.interval_list \
     -R GRCh37.fasta \
     --bqsr-recal-file subsample.recal.table
```

```bash
mv subsample.recal.bai subsample.recal.bam.bai
```

```bash
samtools stats subsample.recal.bam > subsample.recal.stats
```

### 5. Germline Variant Calling

### Haplotype Caller
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

### GenotypeVCFs
```bash
gatk --java-options -Xmx2g \
     GenotypeGVCFs \
     -R GRCh37.fasta \
     -L ../assets/exome.bed.interval_list \
     -D ../assets/dbsnp_138.b37.vcf.gz \
     -V subsample.g.vcf \
     -O subsample.vcf
```

### 6. Select SNPS/INDELS
```bash
gatk SelectVariants \
     -R GRCh37.fasta \
     -V subsample.vcf \
     -O subsample.snps.vcf.gz \
     -select-type SNP
```

```bash
gatk SelectVariants \
     -R GRCh37.fasta \
     -V subsample.vcf \
     -O subsample.indels.vcf.gz \
     -select-type INDEL
```

### 7. Filter SNPS
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

### 8. Filter INDELS
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

### 9. Merge Filtered VCFs
```bash
gatk MergeVcfs \
     -I subsample.filt_indels.vcf \
     -I subsample.filt_snps.vcf \
     -O filtered_sample.vcf
```

### 10. Annotate Variants
Run in your own time, this step takes a long time. Usually we submit this job to SLURM with much more resources than we have initially requested.

```bash
snpEff GRCh37.75 \
       -csvStats subsample.snpEff.csv \
       -nodownload \
       -dataDir /data/snpEff/ \
       -canon \
       -v filtered_sample.vcf > subsample.snpEff.ann.vcf
```
