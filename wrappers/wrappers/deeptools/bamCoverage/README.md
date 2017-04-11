# Deeptools bamCoverage

bamCoverage is a tool to produce bigWigs or bedGraphs. See
[website](http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html#usage-examples-for-rna-seq)
for details.

## Examples

Minimal bigwig usage.

```python
rule bamCoverage:
    input:
        bam='sample1.bam',
        bai='sample1.bam.bai'
    output: 'sample1.bw'
    threads: 2
    wrapper:
        'file://path/to/wrapper'
```

Minimal bedgraph usage.

```python
rule bamCoverage:
    input:
        bam='sample1.bam',
        bai='sample1.bam.bai'
    output: 'sample1.bg'
    params:
        extra='--outFileFormat bedgraph'
    threads: 2
    wrapper:
        'file://path/to/wrapper'
```


Common usage: only output reverse strand, filter reads with mapping quality
less than 20, ignore duplicate reads, smooth by taking average across 10 bases.

```python
rule bamCoverage:
    input:
        bam='sample1.bam',
        bai='sample1.bam.bai'
    output: 'sample1_forward.bw'
    params:
        extra='--filterRNAstrand reverse --minMappingQuality 20 --ignoreDuplicates --smoothLength 10'
    threads: 2
    wrapper:
        'file://path/to/wrapper'
```

## Input
`bam`: A sorted bam file.

`bai`: Bam file must be indexed.

## Output
bamCoverage will output a bigWig by default. It can also output a bedgraph is specified (see example).

## Threads
The wrapper specifically sets `--threads` to `{snakemake.threads}`, which defaults to 1 if
not specified.

## Params
`extra`: passed verbatim to `bamCoverage`.
