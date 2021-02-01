# Wrapper for dupRadar

dupRadar provides an easy way to distinguish between artifactual vs natural
duplicate reads in RNA-Seq data. Prior to dupRadar only global duplication rates
were used and they don't take into account the effect of gene expression levels. 
dupRadar relates *duplication rates* and *length normalized read counts* of every
gene to model the dependency of both variables. 

[Link to homepage](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html)

[Link to manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)

## Example

Single-end, not stranded:

```python
rule dupRadar:
    input:
       bam='sample1.bam',
       annotation='dm6.gtf',
    output:
        density_scatter='sample1.density_scatter.png',
        expression_histogram='sample1.expression_histogram.png',
        expression_boxplot='sample1.expression_boxplot.png',
        expression_barplot='sample1.expression_barplot.png',
        multimapping_histogram='sample1.multimapping_histogram.png',
        dataframe='sample1.dupradar.tsv'
    wrapper:
        wrapper_for('dupRadar')
```

Paired-end, stranded:

```python
rule dupRadar:
    input:
       bam='{sample_dir}/{sample}/{sample}.cutadapt.hisat2.unique.sort.dedup.bam',
       annotation='annotations/dm6.gtf',
    output:
        density_scatter='sample1.density_scatter.png',
        expression_histogram='sample1.expression_histogram.png',
        expression_boxplot='sample1.expression_boxplot.png',
        expression_barplot='sample1.expression_barplot.png',
        dataframe='sample1.dupradar.tsv'
    params:
        paired=True,
        stranded=True
    wrapper:
        wrapper_for('dupRadar')
```

## Input
* `bam`: BAM file with mapped reads has to be duplicate marked using either
  Picard or BamUtil

* `annotation`: GTF file contaning features to count the reads falling on the
  features.

## Output
Output plots are described in the [dupRadar
vignette)[http://bioconductor.org/packages/release/bioc/vignettes/dupRadar/inst/doc/dupRadar.html].
See that page for descriptions of outputs and how to interpret them.

* `density_scatter`: expression vs percent duplication
* `expression_boxplot`: expression vs percent duplication, binned into boxes
* `expression_histogram`: standard histogram of expression (RPKM)
* `expression_barplot`: percentage duplication in 5% expression bins.
* `multimapping_histogram`: histogram showing fraction of reads coming from
  multimapping reads
* `dataframe`: results from `analyzeDuprates` saved as a TSV for downstream
  analysis. Following the vignette, we also add the fraction of multimappers in
  each gene as the column `mhRate`.
* `model`: Slope and intercept of the dupsExpFit
* `curve`: Simplified curve of the GLM for downstream plotting

## Threads
Threads are passed to dupRadar and are in turn passed to featureCounts, which
it calls automatically.

## Params
* `paired`: True | False. Default False.
* `stranded`: True | False | "reverse". Default False.
