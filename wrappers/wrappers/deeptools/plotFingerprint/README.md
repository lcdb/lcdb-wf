# Deeptools plotFingerprint

Evaluates enrichment of a ChIP-seq experiment.

http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html

## Examples

When no `control:` is provided in the input, all BAMs are provided as the
`--bamFiles` argument. Note that the filetype of the plot is determined by the
extension and handled as appropriately by `plotFingerprint`.


```python
rule fingerprint:
    input:
        bams=['ip1.bam', 'ip2.bam', 'input1.bam']
    output:
        plot='fingerprint.png'
    threads: 2
    wrapper:
        'file://path/to/wrapper'
```

Provde a BAM for the `control` in order to have `plotFingerprint` calculate the
metric file. NOTE: if `metrics` is provided, the `--outQualityMetrics` argument
will be automaticaly added so you don't need to include it in the `extra`
params. However, if you do, the `extra` params will take priority.

Note that if you want the output metrics, you probably want to specify the
`control` named input as shown here. This will be added to the call as
`--JSDsample {snakemake.input.control}`

```python
rule fingerprint:
    input:
        bams=['ip1.bam', 'ip2.bam'],
        control='input1.bam'
    output:
        plot='fingerprint.png',
        metrics='fingerprint-metrics.tsv',
        raw_counts='fingerprint-counts.tsv'
    threads: 2
    params:
        extra: '--extendReads=300 --skipZeros'
    wrapper:
        'file://path/to/wrapper'
```


## Input
`bams`: A single BAM or a list of BAMs. Must be indexed; you can provide the
bai files as an additional input argument in order to trigger their creation,
but they are not explicitly handled by the wrapper.

`control`: optional single BAM file that `plotFingerprint` will use to
calculate various metrics.

## Output
`plot`: format determined by extension

`metrics`: optional output, adds `--outQualityMetrics={snakemake.output.metrics}` to extras.

raw_outpu


## Threads

The wrapper specifically sets `--threads` to `{snakemake.threads}`, which
defaults [[to]] 1 if not specified.

## Params
`extra`: passed verbatim to `plotFingerprint`.
