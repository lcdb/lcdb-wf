# Wrapper for MultiQC

MultiQC aggregates output from a wide range of bioinformatics tools into
a single interactive report.

Note that MultiQC has expectations about output filenames but this can be
modified using the `--config` argument in `params.extra`. For example, by
default it looks for FastQC ZIP files that end in `_fastqc.zip`, so you might
want to adjust upstream rules' output to match expectations of MultiQC.

## Example

```
rule multiqc
    input: expand('samples/{sample}{suffix}', sample=samples, suffix=['_fastqc.zip', '.cutadapt.log'])
    output:
        html='multiqc_report.html'
    params:
        analysis_directory='samples'
        extra="--config=multiqc_config.yaml"
    wrapper:
        "file://path/to/multiqc"
```
## Input

Technically MultiQC doesn't need any input, but in order for the rule to
trigger when appropriate you'll want to include as input all of the files
you're expecting it to find.

## Output
- `html`: the output HTML file

Note that the wrapper will automatically set the `--outdir` option to be the
dirname of output.html.

## Threads
Threads not supported.

## Params
`analysis_directory` is required. It points to the top-level dir in which
MultiQC should search.

The `--force` argument is always automatically specified, to allow Snakemake to
decide when to overwrite an existing file when dependent files have been
updated.

Additional parameters can be passed to MultiQC verbatim by supplying a string
in params.extra.

## Log
If snakemake.log is supplied, a log will be written there; otherwise to stdout.

