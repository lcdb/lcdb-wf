# SICER

Wraps the `sicer` program to call ChIP-seq peaks on input BED files.

## Examples

Minimal usage. SICER is the best operating piece of hot garbage you'll ever find.
It has a completely fixed set of input parameters it requires, hard-coded genome
data in SICER/lib/GenomeData.py (submit bug report in bioconda if you need
additions), and it can't be run from the same directory at the same time due to
hard coded output filenames. It's a proper mess boss.

```python
rule sicer:
    input:
        ip='ip.bed',
        control='input.bed',
	redundancy_threshold=1,
	window_size=200,
	fragment_size=150,
	effective_genome_fraction=0.75,
	gap_size=600,
	fdr=0.01
    output:
        bed='out/peaks.bed'
    wrapper:
        'file://path/to/wrapper'
```


## Input

`ip`: single BED for IP

`control`: single BED for input

`redundancy_threshold`: cutoff count above which duplicates are removed

`window_size`: SICER resolution; 200 recommended for histones

`fragment_size`: twice the shift from the beginning to the center of a read

`effective_genome_fraction`: percentage of mappable genome; only set it here if you want to override the genome build in config.yaml

`gap_size`: nonnegative integer multiple of window size. used to merge contiguous regions (higher means more liberal merging).

`fdr`: FDR cutoff for calling significant regions.

## Output

`bed`: BED file of called peaks. This is a delicately processed version of `*island.bed` from SICER.

Other files are created, these can be added as additional named outputs for use
by downstream rules, however the wrapper only pays attention to
`snakemake.output.bed`.


## Params
Do not use `extra` for this rule.
