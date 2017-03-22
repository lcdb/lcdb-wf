# Wrapper for atropos
[Atropos](https://atropos.readthedocs.io/en/latest/index.html) is a fork of
[Cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html) which finds and
removes adapter sequences, primers, poly-A tails and other types of unwanted
sequence from your high-throughput sequencing reads.

# Examples

Minimal usage:

```
rule atropos:
    input: fastq='{sample}.fastq'
    output: fastq='{sample}.trim.fastq'
    threads: 4
    wrapper:
        "file://path/to/atropos"
```

Use an adapters file and quality-trim reads to Q20:

```
rule atropos:
    input: fastq='{sample}.fastq'
    output: fastq='{sample}.trim.fastq'
    params: extra="-a file:adapters.fa -q 20"
    threads: 4
    wrapper:
        "file://path/to/atropos"
```

Optionally provide the adapters file as input in order to trigger a re-run if
it has changed. The wrapper only pays attention to `input.fastq`, so adding
another key doesn't affect the wrapper:

```
rule atropos:
    input:
        fastq='{sample}.fastq',
        adapters='adapters.fa'
    output: fastq='{sample}.trim.fastq'
    params: extra="-a file:adapters.fa -q 20"
    threads: 4
    wrapper:
        "file://path/to/atropos"
```

Example of how to use with other output files. Since the wrapper only pays
attention to `output.fastq`, so other output files can be indicated but their
filenames have to be indicated in `params.`:

```
rule atropos:
    input:
        fastq='{sample}.fastq',
        adapters='adapters.fa'
    output:
        fastq='{sample}.trim.fastq',
        short='{sample}.trim.too-short.fastq',
        untrimmed='{sample}.untrimmed.fastq',
    params:
        extra=(
            "-a file:adapters.fa "
            "-q 20 "
            "--too-short-output={sample}.trim.too-short.fastq "
            "--untrimmed-output={sample}.untrimmed.fastq"
        )
    threads: 4
    wrapper:
        "file://path/to/atropos"
```

You can also run in pair-end mode.

```
rule atropos:
    input:
        R1='{sample}_r1.fastq',
        R2='{sample}_r2.fastq',
        adapters='adapters.fa'
    output:
        R1='{sample}_r1.trim.fastq',
        R1='{sample}_r2.trim.fastq'
    params: extra="-a file:adapters.fa -A file:adapters.fa -q 20"
    threads: 4
    wrapper:
        "file://path/to/atropos"
```


## Input

All inputs are FASTQ files, and they can be optionally gzipped.

### Single-end mode:

fastq : single-end FASTQ file

### Paired-end mode:

R1 : Read 1 FASTQ
R2 : Read 2 FASTQ

See examples below for other input options including adapters.

## Output
q
### Single-end mode:

fastq : Trimmed FASTQ file.

### Paired-end mode:

R1 : trimmed R1 FASTQ file
R2 : trimmed R2 FASTQ file

See examples below for other output options.

## Log
If a log file is specified, stdout and stderr will be captured there.

## Threads
One imporvement of atropos over cutadapt is the ability to use treads which are
passed to the `-T` option.

## Params
Additional parameters can be passed to atropos verbatim by supplying a string
in `params.extra`.
