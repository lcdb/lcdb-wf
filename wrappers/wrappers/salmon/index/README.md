# Salmon index

## Examples

Minimal usage:

```python
rule salmon:
    input:
        fasta='transcriptome.fa'
    output: 'salmon_index/hash.bin'
    wrapper:
        'file://path/to/wrapper'
```

## Input
`fasta`: Transcriptome FASTA file.

## Output
One or more files created by `salmon index`.

A salmon index is a directory of files. Here we specify one of the index files
that's created, `hash.bin`. The wrapper will provide the dirname of the file to
salmon. You can provide any of the files created by `salmon index` here as
a string, list, or dict. If you want the indexing rule to trigger from
a downstream `salmon quant` rule, be sure to use the same filename for that
rule's input.

Other files created by `salmon index` that could also be used:

- `hash.bin`
- `txpInfo.bin`
- `sa.bin`
- `refInfo.json`
- `header.json`
- `versionInfo.json`
- `indexing.log`
- `rsd.bin`
- `quasi_index.log`

## Threads
Does not use threads

## Params
`extra`: passed verbatim to `salmon index`
