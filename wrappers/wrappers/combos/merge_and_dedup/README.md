# Merge and deduplicate

Merges BAM files and then deduplicates the output. However if only one BAM file
is created, the file is simply symlinked.

This wrapper is often needed in ChIP-seq to merge technical replicates. The
same fragment could have been sequenced in multiple tech reps, resulting in
duplicate reads in the merged output even though each individual BAM already
had duplicates removed.

This method has an advantage over merging first and then deduping in separate
rules when we want to retain both individual (per tech rep) deduped BAMs as
well as merged deduped BAMs. Since the deduping has already happened once for
each tech rep, we want to avoid doing so again if no merging happens.

## Examples

Minimal usage:

```python
rule merge_and_dedup:
    input: 'a1.bam', 'a2.bam'
    output:
        bam='a-merged.bam',
        metrics='a-merged.bam.metrics'
    wrapper:
        'file://path/to/wrapper'
```

In the following case, a symlink will be created since no merging needs to be
performed on a single file:

```python
rule merge_and_dedup:
    input: 'a1.bam'
    output:
        bam='a-merged.bam',
        metrics='a-merged.bam.metrics'
    wrapper:
        'file://path/to/wrapper'
```


## Input

Single BAM or list of BAMs.

## Output

- `bam`: output bam file
- `metrics`: optional output metrics file. Default is to use
  `{snakemake.output.bam}.metrics`.

## Threads

Threads are passed to `samtools merge`.

## Params

- `samtools_merge_extra`: addtional args passed verbatim to `samtools merge`

- `markduplicates_extra`: addtional args passed verbatim to `markduplicates_extra`

- `java_args`: passed to MarkDuplicates, often used to provide more memory
  (e.g., `-Xmx32g`). Be sure to increase the corresponding rule's memory
  resource to account for the additional allocation
