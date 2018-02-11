.. _rnaseq:

RNA-seq workflow
================

This workflow is used for RNA-seq and RNA-seq-like analysis (like euRNA-seq,
RIP-seq or small RNA-seq).


This workflow can use references created by the references workflow with no
need to run the references separately. This workflow includes:

    - Build a HISAT2 index
    - Build a salmon transcriptome index
    - Download a GTF annotation
    - Convert the GTF to refflat format
    - Trim reads with cutadapt
    - Align with HISAT2
    - Run FastQC on raw, trimmed, and aligned reads
    - Align reads to rRNA using bowtie2 to evaluate rRNA contamination
    - Count reads in genes with featureCounts
    - Run dupRadar and preseq to assess library complexity
    - Check for evidence of cross-contamination using fastq_screen on multiple
      configured genomes
    - Assess transcript coverage with Picard CollectRnaSeqMetrics
    - Build bigWigs (optionally strand-specific) created from BAM files
    - Aggregate QC results using MultiQC. Includes custom tables for library
      sizes and rRNA contamination
    - Run various QC and differential expression. This is performed in an
      RMarkdown file that runs a standard DESeq2 differential expression
      analysis along with diagnostic plots, exported tables of differentially
      expressed genes for each contrast, and downstream GO analysis using
      clusterProfiler. This file is run and rendered into an output HTML file.
    - Construct and upload a track hub of scaled coverage bigWigs for each
      sample that can be viewed in UCSC Genome Browser

The DAG of jobs looks like this:

.. image:: rnaseq.png

Preliminary setup: AnnotationHub
-------------------------------

In order to support multiple genomes, we use the AnnotationHub package from
Bioconductor as a uniform API to access annotation information for GO analysis
and mapping gene ID to gene name and symbol. This uniformity and flexibility is
offset by the need to do some preliminary setup to work around 1) running the
analysis on nodes with no internet acess and 2) AnnotationHub caching issues.

To run the analysis on a node with no internet access, we need to pre-download
the annotation information. It's best to do this in a directory-specific
manner, since Bioconductor updates their annotations frequently. You don't want
to use the default cache location of ``~/.AnnotationHub``, because that could
collide with other experiments that you don't want to update yet.

The solution is to:
1. activate the environment on a machine with internet access
2. change to the top level directory of the repo
3. start R
4. run the following code:

.. code-block:: r

    library(AnnotationHub)
    ah <- AnnotationHub(cache='include/AnnotationHubCache')
    idx <- grepl('Homo sapiens', ah$species) & ah$rdataclass == 'OrgDb'

    # This may give multiple results. In order to choose, you can view more
    # information on the options:
    mcols(ah[idx])

    # Choose the annotation key you want to use. Accessing the OrgDb object
    # triggers a download, and the download will go to the cache directory
    # specified above.
    orgdb <- ah[[annotation_key]]

Parameters and settings
-----------------------
See :ref:`config` for configuring in general; here we describe the
RNA-seq-specific configuration.

:``config/sampletable.tsv``:
    This TSV is loaded into the RMarkdown file, so any information you add here
    will be available in R and can therefore be used in linear models with
    DESeq2 or can be used to otherwise annotate samples.


:``downstream/rnaseq.Rmd``:
    This file runs the differential expression analysis and other downstream
    work. Search for the string ``Assumption``. This indicates places where
    assumptions have been made about the experiment (what columns to use as
    factors, how to group samples, etc).

:``rnaseq_trackhub.py``:
    This script builds the track hub for RNA-seq. Look for the string
    ``ASSUMPTION`` and follow the instructions there for how to customize for
    your particular experiment.


