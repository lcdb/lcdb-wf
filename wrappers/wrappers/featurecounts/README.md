# Wrapper for featureCounts

Read summarizaton/quantification is required for a number of downstream
analyses such as gene expression and histone modification.  The output
is a counts table where the number of reads assigned to each feature in
each library is recorded. *FeatureCounts* can be used to perform this task
much faster than HTSeq without compromising on accuracy.

[Link to homepage](http://bioinf.wehi.edu.au/featureCounts/)

[Link to manual](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf)

## Examples

Minimal example:

```python
rule featurecounts:
    input:
        annotation='annotation.gtf',
        bam='{sample}.bam'
    output: '{sample}.counts'
    wrapper: 'file://wrapper'
```

Using multiple BAMs, and add extra args (stranded counting):

```python
rule featurecounts:
    input:
        annotation='annotation.gtf',
        bam=expand('{sample}.bam', sample=SAMPLES)
    output: '{sample}.counts'
    params: extra='-s 1'
    wrapper: 'file://wrapper'
```

Adding the summary files to output so that they enter the DAG:

```python
rule featurecounts:
    input:
        annotation='annotation.gtf',
        bam=expand('{sample}.bam', sample=SAMPLES)
    output:
        counts='{sample}.counts',
        summary='{sample}.counts.summary'
    wrapper: 'file://wrapper'
```

## Input
* `bam`: One or more BAM or SAM files
* `annotation`: GTF or GFF file

## Output
* `counts`: The full output of featureCounts, containing all feature related
information such as meta-feature (gene) names, sub-features (e.g., exon)
coordinates, and counts.

`featureCounts` also outputs a summary file containing the number of reads in
various classes (assigned, ambiguous, unmapped, etc). It is automatically saved
alongside the `counts` file, with a `.summary` suffix. If you need to use it
for other downstream rules, it can be added to the output files but will be
ignored by the wrapper (see example above).


## Threads
Threads will be passed to the `-T` parameter.

## Params
* `extra`: a string that will be passes verbatim to featureCounts.
