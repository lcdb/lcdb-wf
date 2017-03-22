# Wrapper for hisat2-build
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) is a fast and
sensitive alignment program for mapping next-generation sequencing reads (both
DNA and RNA) to a population of human genomes (as well as to a single reference
genome).

## Examples

Minimal usage:

```
from lcdblib.snakemake import aligners
rule hisat2_build:
    input: fasta='{sample}.fastq'
    output: index=aligners.hisat2_index_from_prefix('/path/to/index')
    wrapper:
        "file://path/to/wrapper"
```
## Input

fasta: FASTA genome sequence

## Output

index: hisat2 index files.

HISAT2 index files have the format `prefix.N.ht2` where N ranges from 1 to 8.
`hisat2-build` only takes the prefix as an input argument, but snakemake works
with files. So the wrapper figures out what the prefix should be based on the
provided output index files. `output.index` can therefore be just one of the
index files, or a list of them.

The index files can be created by, e.g.,

```python
from lcdblib.snakemake import aligners
indexes = aligners.hisat2_index_from_prefix('/path/to/index')
```

or:

```python
indexes = [
    '{prefix}.{n}.ht2'.format(prefix='/path/to/index', n=n)
    for n in range(1, 9)
]
```

## Log
If a log file is specified, stdout and stderr will be captured there.

## Threads
Ignores threads

## Params
Additional parameters can be passed verbatim by supplying a string in
`params.extra`.


