# Wrapper for Picard MarkDuplicates

## Examples

Minimal usage:

```python
rule markduplicates:
    input:
        bam='{sample}.bam',
    output:
        bam='{sample}.dupsmarked.bam',
        metrics='{sample}.dupmetrics'
    wrapper: 'file://wrapper'
```

Allocate 64GB of memory for Java and remove duplicates instead of just marking them:

```python
rule markduplicates:
    input:
        bam='{sample}.bam',
    output:
        bam='{sample}.dupsmarked.bam',
        metrics='{sample}.dupmetrics'
    wrapper: 'file://wrapper'
    params:
        java_args="-Xmx64g,
        extra: "REMOVE_DUPLICATES=true"
    wrapper: 'file://wrapper'
```

## Input

* `bam`: BAM file

## Output

- `metrics`: a metrics file containing general metrics


## Threads
Threads not supported.

## Params
* `extra`: passed verbatim.
* `java_args` can be used to set Java memory, e.g., `java_args="-Xmx64g"`

## Example output

Metrics file example looks like this:

```
## htsjdk.samtools.metrics.StringHeader
# picard.sam.markduplicates.MarkDuplicates INPUT=[sample1.bam] OUTPUT=sample1.dupsmarked.bam METRICS_FILE=sample1.dupmetrics.txt    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag REMOVE_DUPLICATES=false ASSUME_SORTED=false DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json
## htsjdk.samtools.metrics.StringHeader
# Started on: Thu Dec 29 22:06:14 EST 2016

## METRICS CLASS        picard.sam.DuplicationMetrics
LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED     SECONDARY_OR_SUPPLEMENTARY_RDS  UNMAPPED_READS  UNPAIRED_READ_DUPLICATES        READ_PAIR_DUPLICATES    READ_PAIR_OPTICAL_DUPLICATES    PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
Unknown Library 2166    0       54      242     100     0       0       0.046168
```
