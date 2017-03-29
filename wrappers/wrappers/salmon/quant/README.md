# Salmon quant

## Examples

Minimal usage, single-end. Note the `--libType A` in params:extra.

```python
rule salmon_quant:
    input:
        index='salmon_index/hash.bin',
        unmatedReads='sample1.fastq.gz'
    output:
        'sample1/salmon/quant.sf'
    params:
        extra='--libType A'
    wrapper:
        'file://path/to/wrapper'
```
Minimal usage, paired-end:

```python
rule salmon_quant:
    input:
        index='salmon_index/hash.bin',
        mates1='sample1_R1.fastq.gz',
        mates2='sample1_R2.fastq.gz',
    output:
        'sample1/salmon/quant.sf'
    params:
        extra='--libType A'
    wrapper:
        'file://path/to/wrapper'
```

## Input

`index`: One or more files created by `salmon index`. The wrapper will use the
dirname of these as the index dir.

`unmatedReads`: FASTQ file. Only used for single-end reads

`mates1`, `mates2`: FASTQ files. Only used for paried-end reads.

## Output
Any files created by `salmon quant`. Salmon creates an output directory with at
least `quant.sf` in it; the wrapper will use the dirname of the output file[s]
as the `--output` argument to salmon. There are other possibilities depending
on the arguments -- see the [output format
section](http://salmon.readthedocs.io/en/latest/file_formats.html#fileformats)
of the Salmon docs for more info.

## Threads
Salmon by default will detect the number of available threads. The wrapper
specifically sets `--threads` to `{snakemake.threads}`, which defaults to 1 if
not specified.

## Params
`extra`: passed verbatim to `salmon index`. You'll probably want `--libType A`
   which auto-detects the library type.

