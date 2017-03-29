# Salmon

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
Transcriptome FASTA file.

## Output
A salmon index is a directory of files. Here we specify one of the index files
that's created, `hash.bin`. The wrapper will provide the dirname of the file to
salmon. You can provide any of the files created by `salmon index` here, as
long as you use at least one of the same filenames in a downstream `salmon
quant` rule so it will trigger the index building.

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
