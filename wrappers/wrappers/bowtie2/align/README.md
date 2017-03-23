# Wrapper for bowtie2

This wrapper has 3 stages:

- align with bowtie2 (output is temporary SAM file)
- convert SAM to BAM (output is BAM)
- sort the BAM (output is sorted BAM)

The `params` section of a rule can pass extra arguments to each stage via the
`*_extra` values. See the *Params* section below.

## Examples

This example runs bowtie2 on single-end reads in very sensitive mode:

```python
from lcdblib.snakemake import aligners
rule bowtie2:
    input:
        fastq="mapped/{sample}_R1.fastq.gz",
        index=aligners.bowtie2_index_from_prefix('/data/assembly/assembly')
    output:
        bam="mapped/{sample}.sorted.bam"
    params: bowtie2_extra="--very-sensitive"
    wrapper:
        "file://wrapper"
```


Example of paired-end mode:

```python
from lcdblib.snakemake import aligners
rule bowtie2:
    input:
        R1="{sample}_R1.fastq.gz",
        R2="{sample}_R2.fastq.gz",
        index=aligners.bowtie2_index_from_prefix('/data/assembly/assembly')
    output:
        bam="mapped/{sample}.sorted.bam"
    wrapper:
        "file://wrapper"
```


When aligning, use 8 threads.  When creating the BAM, filter out unmapped
reads. When sorting, use 4 threads and use the highest compression level.

```python
from lcdblib.snakemake import aligners
rule bowtie2:
    input:
        fastq="mapped/{sample}_R1.fastq.gz",
        index=aligners.bowtie2_index_from_prefix('/data/assembly/assembly')
    output:
        bam="mapped/{sample}.sorted.bam"
    threads: 8
    params:
        samtools_view_extra="-F 0x04",
        samtools_sort_extra="--threads 4 -l 9"
    wrapper:
        "file://wrapper"
```

Map paired-end reads, but return a name-sorted BAM rather than
a coordinate-sorted BAM.

```python
from lcdblib.snakemake import aligners
rule bowtie2:
    input:
        R1="{sample}_R1.fastq.gz",
        R2="{sample}_R2.fastq.gz",
        index=aligners.bowtie2_index_from_prefix('/data/assembly/assembly')
    output:
        bam="mapped/{sample}.sorted.bam"
    params:
        samtools_sort_extra='-n'
    wrapper:
        "file://wrapper"
```
## Input

### Single-end mode

`fastq`: a single FASTQ file.  Can be gzipped.

`index`: List of the `*.bt2` index files. See the examples for a good way to
create this list. The wrapper will figure out the prefix to provide to bowtie2
based on these files.


### Paired-end mode

`R1`: read 1 FASTQ file, can be gzipped
`R2`: read 2 FASTQ file, can be gzipped

`index`: List of the `*.bt2` index files. See the examples for a good way to
create this list. The wrapper will figure out the prefix to provide to bowtie2
based on these files.

## Output

`bam`: Coordinate-sorted BAM file of aligned reads

## Threads

Threads are passed to bowtie2 as the `--threads` parameter. `samtools sort` also
accepts a `--threads` parameter, but since it uses 768MB memory per thread,
this is set manually using the `samtools_extra` param described below

## Params

`bowtie2_extra` can be a string of arbitrary parameters that will be passed verbatim to
bowtie2.

`samtools_view_extra` can be a string of arbitrary parameters passed verbatim
to `samtools view`. Note that the wrapper already uses `-Sb` to convert SAM to
BAM.

`samtools_sort_extra` can be a string of parameters passed verbatim to
`samtools sort`. Note that the wrapper already uses `-O=BAM -o={output.bam}`.

