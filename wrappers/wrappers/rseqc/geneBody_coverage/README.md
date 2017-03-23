# rseqc geneBody Coverage

[website](http://rseqc.sourceforge.net/#genebody-coverage-py)

## Examples

### Minimal usage:

```python
rule geneBody_coverage:
    input: bam = 'sample1.sort.bam',
           bai = 'sample1.sort.bam.bai',
           bed = 'dm6.bed'
    output: txt='sample1.geneBodyCoverage.txt',
            r='sample1.geneBodyCoverage.r',
            img='sample1.geneBodyCoverage.pdf'
    log: 'sample1.geneBodyCoverage.log'
    wrapper:
        'file://path/to/wrapper'
```

### Additional Arguments:

```python
rule geneBody_coverage:
    input: bam = 'sample1.sort.bam',
           bai = 'sample1.sort.bam.bai',
           bed = 'dm6.bed'
    output: txt='sample1.geneBodyCoverage.txt',
            r='sample1.geneBodyCoverage.r',
            img='sample1.geneBodyCoverage.png'
    params: 
        extra = '-f png'
    log: 'sample1.geneBodyCoverage.log'
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
`txt`:
    Results table.
`r`:
    R script used to generate figure.
`img`:
    image output version of the figure [pdf is default] to change specify `-f`
    option in `params.extra`

## Threads
Does not use threads

## Params
`extra`:
    Command line parameters that will be passed directly to rseqc
    geneBody-coverage.py See website for possible options.

