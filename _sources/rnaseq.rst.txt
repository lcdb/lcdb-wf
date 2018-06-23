.. _rnaseq:

RNA-seq workflow
================

This workflow is used for RNA-seq and RNA-seq-like analysis (like euRNA-seq,
RIP-seq or small RNA-seq).

This workflow can use references created by the references workflow with no
need to run the references workflow separately. This workflow performs the
following tasks:

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
    - Optionally merge bigWigs as defined by config
    - Aggregate QC results using MultiQC. Includes custom tables for library
      sizes and rRNA contamination
    - Run various QC and differential expression. This is  performed in an
      RMarkdown file that runs a standard DESeq2 differential expression
      analysis along with diagnostic plots, exported tables of differentially
      expressed genes for each contrast, and downstream GO analysis using
      clusterProfiler. This file is run and rendered into an output HTML file.
    - Construct and upload a track hub of scaled coverage bigWigs for each
      sample that can be viewed in UCSC Genome Browser

The DAG of jobs looks like this:

.. image:: rnaseq.png

To run the workflow on your own data:

- Install the dependencies (:ref:`getting-started`)
- Optionally run the workflow on provided test data (:ref:`running-the-tests`)
- Search the Snakefile for the string ``NOTE`` and edit rules as described by those comments
- Edit the sampletable to reflect your sample IDs and their FASTQ files (:ref:`sampletable`)
- Edit the config file to reflect the references you want to use (:ref:`references-config`)
- Edit the config file to reflect the optional bigWig merging you would like (:ref:`cfg-merged-bigwigs`)

If you are running on a local machine, from the ``workflows/rnaseq`` directory
and with the environment activated, this will do a dry run:

.. code-block:: bash

    snakemake -n

If no errors are displayed, run the workflow (assuming 8 CPUs):

.. code-block:: bash

    snakemake -j 8 --use-conda

If you are running on a cluster other than NIH's Biowulf, edit the settings in
the ``workflows/rnaseq/config/clusterconfig.yaml`` to reflect the settings for
your cluster.

Edit ``include/WRAPPER_SLURM`` to reflect the environment name that you used.

From the ``workflows/rnaseq`` directory:

.. code-block:: bash

    sbatch ../../include/WRAPPER_SLURM

