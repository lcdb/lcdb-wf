# Wrapper for bowtie2-build
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an
ultrafast and memory-efficient tool for aligning sequencing reads to long
reference sequences.

## Examples

Minimal usage:

```
from lcdblib.snakemake import aligners
rule bowtie2_build:
    input: fasta='{sample}.fastq'
    output: index=aligners.bowtie2_index_from_prefix('/path/to/index')
    wrapper:
        "file://path/to/wrapper"
```
## Input

fasta: FASTA genome sequence

## Output

index: bowtie2 index files.

bowtie2 index files have the format `prefix.N.bt2` where N ranges from 1 to 4.
`bowtie2-build` only takes the prefix as an input argument, but snakemake works
with files. So the wrapper figures out what the prefix should be based on the
provided output index files. `output.index` can therefore be just one of the
index files, or a list of them.

The index files can be created by, e.g.,

```python
from lcdblib.snakemake import aligners
indexes = aligners.bowtie2_index_from_prefix('/path/to/index')
```

or:

```python
indexes = [
    '{prefix}.{n}.bt2'.format(prefix='/path/to/index', n=n)
    for n in range(1, 5)
]
```

## Log
If a log file is specified, stdout and stderr will be captured there.

## Threads
Ignores threads

## Params
Additional parameters can be passed verbatim by supplying a string in
`params.extra`.


