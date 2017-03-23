# rseqc Infer Experiment

[website](http://rseqc.sourceforge.net/#infer-experiment-py)

## Examples

### Minimal usage:

```python
rule infer_experiment:
    input: bam = 'sample1.bam',
           bed = 'dm6.bed'
    output: 
        txt = 'sample1.infer_experiment.txt'
    log: 'sample1.infer_experiment.log'
    wrapper:
        'file://path/to/wrapper'
```

### Extra Arguments:

```python
rule infer_experiment:
    input: bam = 'sample1.bam',
           bed = 'dm6.bed'
    output: 
        txt = 'sample1.infer_experiment.txt'
    params:
        extra = '-q 20'
    log: 'sample1.infer_experiment.log'
    wrapper:
        'file://path/to/wrapper'
```

## Input
`bam`:
    Input alignment file in BAM format.

`bed`:
    Reference gene model in bed format.

## Output
`txt`:
    Contains infer_experiment output..

## Threads
Does not use threads

## Params
`extra`:
    Command line parameters that will be passed directly to rseqc infer_experiment.py.
    See website for possible options.
