# MA5112 Group 1 Feedback

Totally understandable if you never want to look at nextflow again, but if you do... here are some suggestions for your code.

***

```nextflow
params.input = ''
if(!params.input) exit 1, "error: No input data provided."
Channel.fromPath(params.input)
       .into{ reads_fastqc_ch; reads_trimming_ch }


params.genome = ''
if(!params.genome) exit 1, "error: No reference genome provided."
Channel.fromPath(params.genome)
       .set{ reference_ch }

outdir = "."
```

Very small issue, but `outdir` is a hardcoded variable. Alternatively, you should let users define this via the command line i.e `params.outdir` in the code, then `--outdir "some_path/"` when running on cmd line.

I never covered this in the tutorials, but to declutter your script it is common practice to put your parameters in a configuration profile. So in your `nexflow.config` file, you might use:

```nextflow
params{
  outdir = null
  genome = null
  input = null
}
```

***

Nice work on the help message, appreciated that!

***

For `bismark`, it looks like you can toggle between `Bowtie` and `Hisat2` by providing a string when running the command, which is a perfect opportunity to use a parameter to define which aligner to use. Check out the code below: 

```nextflow

params.bismark_aligner = 'Hisat2'

process bismark_genome_preparation{

    publishDir "${outdir}/results_nextflow", mode:'copy'

    input:
    file(reference) from reference_ch

    output:
    file"Bismark_Index" into index_ch, index2_ch

    script:
    """
    mkdir Bismark_Index
    bismark_genome_preparation ${params.bismark_aligner} $reference
    cp -r $reference Bismark_Index/
    """
}
```
