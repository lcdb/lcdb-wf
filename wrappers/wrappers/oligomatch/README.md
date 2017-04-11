# oligomatch


## Examples

Minimal usage:

```python
rule oligomatch:
    input:
        oligos='HindIII.fa'
        sequence='dm6.fa'
    output: 'restriction_sites.bed'
    wrapper:
        'file://path/to/wrapper'
```

## Input

`oligos` and `sequence` can be `.fa`, `.nib`, or `.2bit` format.

## Threads
Does not use threads

## Params
oligoMatch does not take any additional parameters.
