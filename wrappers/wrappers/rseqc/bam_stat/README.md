# rseqc Infer Experiment

[website](http://rseqc.sourceforge.net/#bam-stat-py)

## Examples

### Minimal usage:

```python
rule bam_stat:
    input: 
        bam = 'sample1.bam'
    output: 
        txt='sample1.bam_stat.txt'
    wrapper:
        'file://path/to/wrapper'
```

### Extra Arguments:

```python
rule bam_stat:
    input: 
        bam = 'sample1.bam'
    output: 
        txt='sample1.bam_stat.txt'
    params:
        extra='-q 20'
    wrapper:
        'file://path/to/wrapper'
```

## Input

`bam`:
    Input alignment file in BAM format.

## Output

`txt`: 
    Text file containing basic stats.

## Threads
Does not use threads

## Params
`extra`:
    Command line parameters that will be passed directly to rseqc bam_stat.py.
    See website for possible options.

## Log
Cannot use log, bam_stat.py writes output to STDERR which is captured in the output.txt option.
