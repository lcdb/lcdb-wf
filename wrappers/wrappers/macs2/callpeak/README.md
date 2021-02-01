# MACS2

Wraps the `macs2 callpeak` subprogram to call ChIP-seq peaks on input BAM
files.

## Examples

Minimal usage. MACS2 outputs a whole directory; this directory is the dirname
of `output.bed`. Note the specification of the genome size in `params.extra`.

```python
rule macs2:
    input:
        treatment='ip.bam',
        control='input.bam',
        chromsizes='dm6.chromsizes'
    output:
        bed='out/peaks.bed'
    extra: '-g dm'
    wrapper:
        'file://path/to/wrapper'
```

MACS2 supports multiple ip and input samples (they are concatenated). This also
shows broad peak-calling, asks MACS2 to create scaled bedgraphs, and adds them as
output files so downstream rules can use them:

```python
rule macs2:
    input:
        treatment=['ip1.bam', 'ip2.bam'],
        control=['input1.bam', 'input2.bam'],
        chromsizes='dm6.chromsizes'
    output:
        bed='out/peaks.bed'
    params: extra='-g dm --bdg --SPMR --broad'
    wrapper:
        'file://path/to/wrapper'
```

## Input

`treatment`: single BAM or list of BAMs for IP

`control`: single BAM or list of BAMs for input

`chromsizes`: Chromsizes table, used to ensure peak boundaries do not extend
outside of chromosome limits.

## Output

`bed`: BED file of called peaks. This is symlinked from the
`*_peaks.narrowPeak` or `*_peaks.broadPeak` file created by MACS2.

Other files are created, these can be added as additional named outputs for use
by downstream rules, however the wrapper only pays attention to
`snakemake.output.bed`.


## Params
Additional params in `extra` will be passed verbatim to `macs2 callpeak`.
