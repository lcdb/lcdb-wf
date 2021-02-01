# spp

Wraps the [`spp`](http://compbio.med.harvard.edu/Supplements/ChIP-seq/) peak-caller.

This is a rather complicated wrapper. See input and output sections below for
details.


## Examples

Minimal usage:

```python
rule spp:
    input:
      ip="ip.bam",
      control="control.bam",
      chromsizes='dm6.chromsizes'
    output: "peaks.bed"
    wrapper:
        'file://path/to/wrapper'
```

Specify parameters (see below for options):


```python
rule spp:
    input:
      ip="ip.bam",
      control="control.bam",
      chromsizes='dm6.chromsizes'
    output: "peaks.bed"
    params: block={'fdr': 0.1}

    wrapper:
        'file://path/to/wrapper'
```

Specify additional output files:

```python
rule spp:
    input:
        ip="ip.bam",
        control="control.bam",
        chromsizes='dm6.chromsizes'
    output:
        bed="peaks.bed"
        enrichment_estimates="enrichment_est.bedgraph",
        smoothed_enrichment_mle="enrichment_mle.bedgraph",
        rdata="image.RData"
    params: block={'fdr': 0.1}
    log: "spp.log"
```

The works, with multiple replicate BAMs to be merged, keeping the tempfiles,
increasing the memory available to MarkDuplicates, all the output files,
adjusting spp params, and using 8 threads for merging and duplicates removal:


```python
rule spp:
    input:
        ip=["ip.bam", "ip2.bam"],
        control=["control.bam", "control2.bam", "control3.bam"],
        chromsizes='dm6.chromsizes'
    output:
        bed="peaks.bed"
        enrichment_estimates="enrichment_est.bedgraph",
        smoothed_enrichment_mle="enrichment_mle.bedgraph",
        rdata="image.RData"
    log: 'spp.log'
    threads: 8
    params:
        block={'fdr': 0.1, 'bins': 10},
        java_args='-Xmx64g'
        keep_tempfiles=True
    log: "spp.log"
```

## Input

`ip`, `control`: BAM files. Duplicates should already be removed.

`chromsizes`: Chromsizes table, used to ensure peak boundaries do not extend
outside of chromosome limits.

SPP itself only supports a single BAM file for IP and a single BAM file for
control.  However, to support the common case of pooling replicates to gain
coverage, this wrapper does handle multiple BAMs.

If more than one BAM is provided for either IP or control, the BAMs are merged
and then duplicates are removed from the merged file (to handle reads that
occur  in both replicates, which would otherwise cause spp to complain) are
then removed using MarkDuplicates. This merged, deduped BAM is then provided to
SPP.

The merged BAM, merged-and-deduped BAM, and metrics file (from MarkDuplicates)
are created as temp files. The temp filenames are indicated in the log. If you
need these for debugging, set `params: keep_tempfiles=True` to keep them.

## Output

The only required output is `bed`. Others, if specified, will trigger their
respective creation.

`bed`: narrowPeak format.

`smoothed_enrichment_mle`: BEDGRAPH file (even though SPP calls it a "WIG") of
smoothed enrichment using the `smoothed.enrichment.mle` method from SPP.
Optional, if not specified it will not be created.

`enrichment_estimates`: BEDGRAPH file (even though SPP calls it a "WIG") of
enrichment estimates using the `get.conservative.fold.enrichment.profile`
function from SPP. Optional, if not specified will not be created.

`rdata`: Saves an image of the workspace. Handy for debugging. Optional, if not
specified will not be created.

An R script named after the BED file (`{snakemake.output.bed}.R`), will be
written to the output directory. This can be run from the same directory as the
snakefile was run from for debugging purposes.

## Threads
We do not run SPP in parallel mode due to trouble with running the `snow`
library on clusters (it seems to crash unexpectedly and intermittently).
However, for multiple BAMs, we pass the threads to samtools and MarkDuplicates.

## Params

### wrapper params

`keep_tempfiles`: bool; if True then tempfiles created by merging and deduping
replicate BAMs will be retained for debugging purposes.

`java_args`: str; additional args provided to picard, e.g., `java_args="-Xmx64g"`

### spp params

Since SPP doesn't have a command-line interface, we can't use the "extras="
mechanism to pass params verbatim. Instead, the R script created by the wrapper
supports the following parameters, provided as keys to the `block` param to
make it easier to work with the chipseq config format. For example:

```python
params:
    block={'bins': 5, 'fdr': 0.1},
    java_args='-Xmx64g'
```

`srange`: tuple; controls the range of lags over which to calculate
cross-correlation. Default is `(50, 500)`

`bins`: integer; controls how the binding characteristics will be binned. Default
is `5`.

`tecfilter`: bool; passed to `find.binding.positions` function. Default is True;
set to False to prevent the exclusion of large regions with higher input than
expected.

`remove_anomalies`: bool; enable/disable the remove.tag.anomalies step. Defualt
is False (do not remove anomalies). Setting to True can increase the time
dramatically.

`fdr`: float; false discovery rate when calling peaks. Default is `0.05`.

`whs`: int. window half-size. Used if the auto-calculated
`binding.characteristics` is NA. Default is `500`.

`zthr`: float. Z threshold used when adding broad regions. Default is `3`.

`bandwidth`: int. Bandwith for smoothing WIG file. Default is `200`.

`step`: int; step size for smoothing WIG file. Default is `100`.
