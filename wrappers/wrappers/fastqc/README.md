# Wrapper for FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) performs
quality control for high-throughput sequencing data.

## Input
FASTQ, SAM, or BAM file. FastQC will auto-detect, but you can also use
`--format` and one of bam, sam, bam_mapped, sam_mapped or fastq in the
params.extra field (see example).

## Output
- html: an html file containing the report for the sample
- zip: a zip file containing the images and text file of results

## Threads
Supports threads, passed in as the `--threads` arg

## Params
Additional parameters can be passed to FastQC verbatim by supplying a string in params.extra.

# Example

```
rule fastqc:
    input: 'samples/{sample}.fastq'
    output:
        html='samples/{sample}.fastqc.html',
        zip='samples/{sample}.fastqc.zip'
    params: extra="--contaminants adapters.tsv --format fastq"
    wrapper:
        "file://path/to/fastqc"
```
