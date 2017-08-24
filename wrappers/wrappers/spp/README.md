# spp

Wraps the [`spp`](http://compbio.med.harvard.edu/Supplements/ChIP-seq/) peak-caller.

## Examples

Minimal usage:

```python
rule spp:
    input:
      ip="ip.bam",
      control="control.bam"
    output: "peaks.bed"
    wrapper:
        'file://path/to/wrapper'
```

Specify parameters (see below for options):


```python
rule spp:
    input:
      ip="ip.bam",
      control="control.bam"
    output: "peaks.bed"
    params: block={'fdr': 0.1}
    
    wrapper:
        'file://path/to/wrapper'
```

```python
rule spp:
    input:
        ip="ip.bam",
        control="control.bam"
    output:
        bed="peaks.bed"
        enrichment_estimates="enrichment_est.bedgraph",
        smoothed_enrichment_mle="enrichment_mle.bedgraph",
        rdata="image.RData"
    log: "spp.log"
```

## Input

BAM files. Duplicates should already be removed.

## Output

bed: narrowPeak format.

smoothed_enrichment_mle: BEDGRAPH file (even though SPP calls it a "WIG") of
smoothed enrichment using the `smoothed.enrichment.mle` method from SPP.
Optional, if not specified it will not be created.

enrichment_estimates: BEDGRAPH file (even though SPP calls it a "WIG") of
enrichment estimates using the `get.conservative.fold.enrichment.profile`
function from SPP. Optional, if not specified will not be created.

rdata: Saves an image of the workspace. Handy for debugging. Optional, if not
specified will not be created.

An R script named after the BED file (`{snakemake.output.bed}.R`), will be
written to the output directory. This can be run from the same directory as the
snakefile was run from for debugging.

## Threads
Does not use threads

## Params

Since SPP doesn't have a command-line interface, we can't use the "extras="
mechanism to pass params verbatim. Instead, the R script created by the wrapper
supports the following parameters:

srange: tuple; controls the range of lags over which to calculate cross-correlation. Default is `(50, 500)`

bins: integer; controls how the binding characteristics will be binned. Default is `5`.

tecfilter: bool; passed to `find.binding.positions` function. Default is True;
set to False to prevent the exclusion of large regions with higher input than
expected.

remove_anomalies: bool; enable/disable the remove.tag.anomalies step. Defualt
is False (do not remove anomalies). Setting to True can increase the time
dramatically.

fdr: float; false discovery rate when calling peaks. Default is `0.05`.

whs: int. window half-size. Used if the auto-calculated `binding.characteristics` is NA. Default is `500`.

zthr: float. Z threshold used when adding broad regions. Default is `3`.

bandwidth: int. Bandwith for smoothing WIG file. Default is `200`.

step: int; step size for smoothing WIG file. Default is `100`.
