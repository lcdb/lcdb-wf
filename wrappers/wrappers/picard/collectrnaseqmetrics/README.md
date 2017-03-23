# Wrapper for Picard CollectRnaSeqMetrics



## Examples

**Note: `STRAND` or `STRAND_SPECIFICITY` must be included in `params.extra`.**
From the help for CollectRnaSeqMetrics:

```
STRAND_SPECIFICITY=StrandSpecificity
STRAND=StrandSpecificity
    For strand-specific library prep. For unpaired reads, use
    FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the
    transcription strand.  Required. Possible values: {NONE,
    FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}
```

Minimal usage (unstranded library):

```python
rule collectrnaseqmetrics:
    input:
        bam='{sample}.bam',
        refflat='annotation.refflat'
    output:
        metrics='{sample}.rnametrics'
    params:
        extra="STRAND=NONE"
    wrapper: 'file://wrapper'
```

Allocate 64GB of memory for Java:

```python
rule collectrnaseqmetrics:
    input:
        bam='{sample}.bam',
        refflat='annotation.refflat'
    output:
        metrics='{sample}.rnametrics'
    params:
        extra="STRAND=NONE",
        java_args="-Xmx64g
    wrapper: 'file://wrapper'
```

Output a plot a figure of 5'-3' bias for a single-end stranded library, and
also include the PDF in the output so it enters the DAG:

```python
rule collectrnaseqmetrics:
    input: '{sample}.bam'
    output:
        metrics='{sample}.rnametrics'
        plot='{sample}.bias.pdf'
    params:
        extra='CHART_OUTPUT={sample}.bias.pdf STRAND=FIRST_READ_TRANSCRIPTION_STRAND'
    wrapper: 'file://wrapper'
```
## Input

* `bam`: BAM file
* `refflat`: RefFlat format file. If all you have is a GTF, use UCSC's
gtfToGenePred to convert (`ucsc-gtftogenepred` in bioconda).

## Output

- `metrics`: a metrics file containing general metrics at the top followed by a
histogram of normalized coverage across the average transcript. The output
looks something like this:

```
## htsjdk.samtools.metrics.StringHeader
# picard.analysis.CollectRnaSeqMetrics REF_FLAT=dm6.refflat STRAND_SPECIFICITY=NONE INPUT=sample1.bam OUTPUT=sample1.metrics    MINIMUM_LENGTH=500 RRNA_FRAGMENT_PERCENTAGE=0.8 METRIC_ACCUMULATION_LEVEL=[ALL_READS] ASSUME_SORTED=true STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json
## htsjdk.samtools.metrics.StringHeader
# Started on: Thu Nov 17 10:11:01 EST 2016

## METRICS CLASS>...picard.analysis.RnaSeqMetrics
PF_BASES>...PF_ALIGNED_BASES>...RIBOSOMAL_BASES>CODING_BASES>...UTR_BASES>..INTRONIC_BASES>.INTERGENIC_BASES>...IGNORED_READS>..CORRECT_STRAND_READS>...INCORRECT_STRAND_READS>.PCT_RIBOSOMAL_BASES>PCT_CODING_BASES>...PCT_UTR_BASES>..PCT_INTRONIC_BASES>.PCT_INTERGENIC_BASES>...PCT_MRNA_BASES>.PCT_USABLE_BASES>...PCT_CORRECT_STRAND_READS>...MEDIAN_CV_COVERAGE>.MEDIAN_5PRIME_BIAS>.MEDIAN_3PRIME_BIAS>.MEDIAN_5PRIME_TO_3PRIME_BIAS>...SAMPLE>.LIBRARY>READ_GROUP
86653>..77595>..>...0>..73048>..4370>...177>0>..0>..0>..>...0>..0.941401>...0.056318>...0.002281>...0.941401>...0.842994>...0>..0.710012>...1.263827>...0.984555>...1.259615>...>...>...

## HISTOGRAM>...java.lang.Integer
normalized_position>All_Reads.normalized_coverage
0>..1.388667
1>..1.588762
2>..1.026843
3>..0.997843
4>..0.920145
5>..0.696602
6>..0.642873
7>..0.773903
8>..0.610144


...
```

## Threads
Threads not supported.

## Params
Additional parameters can be passed verbatim. **`STRAND` must be included**,
see above note for details.

* `java_args` can be used to set Java memory, e.g., `java_args="-Xmx64g"
