# Wrapper for Kallisto index

Builds an index for Kallisto.

## Examples

Minimal usage:

```python
rule kallisto_index:
    input: fasta='transcriptome.fa'
    output: index='transcriptome.idx'
    wrapper:
        'file://path/to/wrapper'
```

Use a log and specify where the `stats` file should go:

```python
rule kallisto_index:
    input: fasta='transcriptome.fa'
    output:
        index='transcriptome.idx',
        stats='transcriptome.idx.stats'
    log: 'transcriptome.idx.log'
    wrapper:
        'file://path/to/wrapper'
```

## Input

* `fasta`: transcriptome FASTA file

## Output

* `index`: Filename of kallisto index file.

## Threads
Does not use threads

## Params
`params.extra` are passed verbatim.
