# rseqc Infer Experiment

[website](http://rseqc.sourceforge.net/#tin-py)

## Examples

### Minimal usage:

```python
rule tin:
    input:
        bam='sample1.sort.bam',
        bai='sample1.sort.bam.bai',
        bed='dm6.bed12'
    output: 
        table='sample1.tin.tsv',
        summary='sample1.tin.summary.txt'
    log: 'sample1.tin.log'
    wrapper:
        'file://path/to/wrapper'
```

### Additional Arguments:

```python
rule tin:
    input:
        bam='sample1.sort.bam',
        bai='sample1.sort.bam.bai',
        bed='dm6.bed12'
    output: 
        table='sample1.tin.tsv',
        summary='sample1.tin.summary.txt'
    params:
        extra = '--subtract-background'
    log: 'sample1.tin.log'
    wrapper:
        'file://path/to/wrapper'
```

## Input

`bam`:
    Input alignment file in sorted BAM format.

`bai`:
    Samtools index of BAM file.

`bed`:
    Reference gene model in bed format.

## Output

`table`:
    TSV of results

`summary`:
    Summary of results

## Threads
Does not use threads

## Params
`extra`:
    Command line parameters that will be passed directly to rseqc tin.py.
    See website for possible options.

