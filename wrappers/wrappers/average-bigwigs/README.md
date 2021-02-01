# Average bigWigs

Often we'd like to merge multiple bigWigs together for downstream work
(heatmaps, etc) but there's no single tool to do this. This wrapper runs
`bigWigMerge` on the inputs to sum their values, then uses `awk` to divide by
their values and sort the way bedGraphToBigWig wants them.

The intermediate bedGraph file will be created in ``$TMPDIR``.

## Examples

Minimal usage:

```python
rule average_bigwigs:
    input: 
        bigwigs=[
            'a.bw',
            'b.bw',
            'c.bw'],
        chromsizes='genome.chromsizes'
    output:
        'out.bw'
    wrapper:
        'file://path/to/wrapper'
```

Increase memory used for sorting:

```python
rule average_bigwigs:
    input: 
        bigwigs=[
            'a.bw',
            'b.bw',
            'c.bw'],
        chromsizes='genome.chromsizes'
    output:
        'out.bw'
    params:
        memory='32G'
    wrapper:
        'file://path/to/wrapper'
```

Single bigwig just gets symlinked over.

```python
rule average_bigwigs:
    input: 
        bigwigs='a.bw',
        chromsizes='genome.chromsizes'
    output:
        'out.bw'
    params:
        memory='32G'
    wrapper:
        'file://path/to/wrapper'
```

## Input

List of bigWig files.


## Output

Single bigWig file created by averaging the inputs

## Threads
Does not use threads

## Params

memory: Passed to `sort` as the `-S` argument.
