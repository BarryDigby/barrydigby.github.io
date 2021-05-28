#!/usr/bin/env nextflow

/*
  We used the methylseq main.nf script as a reference for most of this workflow:
  https://github.com/nf-core/methylseq

  Some additional material was also found at /data/bdigby/rugrats when combining channels in step 6 that helped us over the line.
*/

// Input parameters (reference genome and input reads)


// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.input = ''
if(!params.input) exit 1, "error: No input data provided."
Channel.fromPath(params.input)
       .into{ reads_fastqc_ch; reads_trimming_ch }


params.genome = ''
if(!params.genome) exit 1, "error: No reference genome provided."
Channel.fromPath(params.genome)
       .set{ reference_ch }

outdir = "."

def helpMessage() {
    log.info"""
    =============================================
          MA5112 Assignment WGBS
    =============================================
    Dependencies:

    nexflow script
    /data/containers/wgbs/wgbs_final.nf

    reference genome
    /data/containers/wgbs/reference_genome/

    input reads
    /data/containers/wgbs/data/*.fastq

    singularity container
    /data/containers/wgbs.img


    Usage:

    The typical command for running the pipeline:

    nextflow -bg run wgbs_final.nf \\
             --input "/data/containers/wgbs/data/*.fastq" \\
             --genome "/data/containers/wgbs/reference_genome/" \\
             -with-singularity /data/containers/wgbs/wgbs.img

    Arguments:
      --input                           [path] Path to input WGBS data. Use suitable glob pattern to capture all reads and surround in double quotes

      --genome                          [path] Path to directory containing reference genome

      -with-singularity                 [path] Path to singularity container

     """.stripIndent()
}


/*
  Step 0: Preprocessing - Indexing
  Kevin

  Run script only once to prepare the genome of interest for bisulfite alignments.
  Specify a directory containing the genome you want to align your reads against.
  We could choose between Bowtie and HISAT2 aligners and decided on HISAT2, so --hisat2 needs to be included to create a genome index for use with HISAT2.
*/

process bismark_genome_preparation{

    publishDir "${outdir}/results_nextflow", mode:'copy'

    input:
    file(reference) from reference_ch
    // reference genome

    output:
    file"Bismark_Index" into index_ch, index2_ch
    // index_ch used in alignment step
    // index2_ch used in methylation call step

    script:
    """
    mkdir Bismark_Index
    bismark_genome_preparation --hisat2 $reference
    cp -r $reference Bismark_Index/
    """
}

/*

  Step 1: QC on raw data
  Stephanie

  Use Fastqc to QC the raw data. Creates html reports that can be viewed later, though not used in later processes.

  Fastqc Arguments:
  fastqc --quiet (does not display update messages except for errors)

*/

process fastqc{

    publishDir "${outdir}/results_nextflow/fastqc_out", mode:'copy'

    input:
    file reads from reads_fastqc_ch

    output:
    file "*.{html,zip}" into QCresults_ch

    script:
    """
    fastqc --quiet $reads
    """
}

/*

  Step 2: Adapter Trimming
  Hannah

  Use Trim Galore to trim adapters.

  Trim Galore Arguments:
  trim_galore --fastqc

*/

process adapter_trimming{

    publishDir "${outdir}/results_nextflow/trimmed_reads", mode:'copy'

    input:
    file reads from reads_trimming_ch

    output:
    file "*.fq" into trimmed_reads_ch
    file "*.{txt,html,zip}"

    script:
    """
    trim_galore --fastqc $reads
    """
}

/*
   Step 3: Align Reads with Bismark
   Charlie

   Use bismark to perform bisulfite alignment and methylation calling.
   Requires genome of interest and the files to be analysed.

   Bismark Arguments:
      --hisat2
      --genome    (from preprocessing step)
*/


process align_bismark{

    publishDir "${outdir}/results_nextflow/aligned_bismark", mode:'copy'

    input:
    file index from index_ch
    file trimmed_reads from trimmed_reads_ch.collect()

    output:
    file "*.bam" into aligned_reads1_ch, aligned_reads_for_summary_report_ch
    file "*report.txt" into aligned_reads2_ch, aligned_reads_report_for_summary_ch

    script:
    """
    bismark $trimmed_reads --genome "${index}/reference_genome" --hisat2
    """
}


/*
   Step 4: Deduplication
   Luke

   Use deduplicate_bismark to deduplicate the Bismark alignment BAM file and remove all reads but one which align to
   the very same position and in the same orientation.

   deduplicate_bismark Arguments:
                   -s             for single end reads
*/




process deduplicate_bismark{

    publishDir "${outdir}/results_nextflow/deduplicate_output", mode:'copy'

    input:
    file aligned_reads from aligned_reads1_ch.collect()

    output:
    file "*.bam" into deduplicated_reads1_ch
    file "*.deduplication_report.txt" into deduplicated_reads2_ch, deduplicated_reads3_ch


    script:
    """
    deduplicate_bismark $aligned_reads -s
    """
}

/*

  Step 5: Extract Methylation Calls
  Hannah, Kevin
*/
//  Use bismark_methylation_extractor script which takes Bismark result files and extracts the methylation call for every single C analysed

//  bismark_methylation_extractor Arguments:


process bismark_methylation_extractor{

    publishDir "${outdir}/results_nextflow/extract_output", mode:'copy'

    input:
    file deduplicated_reads from deduplicated_reads1_ch.collect()
    file index2 from index2_ch

    output:
    file("*splitting_report.txt") into splitting_report_ch, splitting_summary_ch, splitting_multiqc_ch
    file("*.M-bias.txt") into mbias_report_ch, mbias_summary_ch, mbias_multiqc_ch
    file "*.gz"

    script:
    """
    bismark_methylation_extractor $deduplicated_reads  --comprehensive \\
                --cytosine_report \\
                --bedGraph \\
                --counts \\
                --gzip -s    --genome_folder "${index2}/reference_genome" \
    """
}

/*
  Step 6: Sample Report with Bismark
  Charlie, Luke

  Use bismark2report to generate a html bismark sample report. Optionally further reports of the Bismark suite such as deduplication, methylation extractor (splitting) or M-bias reports can be included.

  bismark2report Arguments:
   --alignment_report
   --dedup_report               optional deduplication report
   --splitting_report           optional splitting report
   --mbias_report               optional M-bias report
*/

// combining multiple channels to be used by bismark2report



aligned_reads_tuple = aligned_reads2_ch
                      .flatten()
                      .map{ it ->  [ it.baseName.substring(0,14), it ] }

dedup_reads_tuple = deduplicated_reads2_ch
                    .flatten()
                    .map{ it -> [ it.baseName.substring(0,14),  it ] }

mbias_tuple = mbias_report_ch
              .flatten()
              .map{ it ->  [ it.baseName.substring(0,14),  it ] }

splitting_report_tuple = splitting_report_ch
                         .flatten()
                         .map{ it -> [ it.baseName.substring(0,14),  it ] }

report_ch = aligned_reads_tuple.join(dedup_reads_tuple).join(mbias_tuple).join(splitting_report_tuple)


process bismark_sample_report{

   publishDir "${outdir}/results_nextflow/sample_report", mode:'copy'

    input:
    tuple val(base), file(hisat2_report), file(dedup_report), file(mbias), file(split_report) from report_ch

    output:
    file "*.{html,txt}" into sample_report_ch

    script:
    """
    bismark2report   --alignment_report $hisat2_report \\
                     --dedup_report $dedup_report \\
                     --splitting_report $split_report \\
                     --mbias_report $mbias

    """
}

/*
  Step 7: Bismark Summary Report
  Kevin, Stephanie

  Combines reports from read alignment, deduplication, and methylation extraction steps for each sample, into one graphical summary html report.
*/
   process bismark_summary_report{

    publishDir "${outdir}/results_nextflow/summary_report", mode:'copy'

    input:
    file ('*') from aligned_reads_for_summary_report_ch.collect()
    file ('*') from aligned_reads_report_for_summary_ch.collect()
    file ('*') from splitting_summary_ch.collect()
    file ('*') from mbias_summary_ch.collect()
    file ('*') from deduplicated_reads3_ch.collect()

    output:
    file "*.html" into bismark_summary_report

    script:
    """
    bismark2summary
    """

}
