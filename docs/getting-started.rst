.. _getting-started:

Getting started
===============
The following steps will install all necessary software and run the example
workflows on a relatively small test data set to ensure that everything runs
correctly. These exampls are run as the automated tests on Travis-CI
(https://travis-ci.org/lcdb/lcdb-wf) to ensure that they are correct.

The example run takes up about 360 MB of space and runs in about 15 mins on
2 cores.

One-time setup
--------------
The following needs to be performed on each system on which you will be running
the workflows.

bioconda (one-time setup)
~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the instructions for setting up `bioconda
<https://bioconda.github.io>`_.  This includes installing `conda` and setting
up the channels in the correct order.

This is required to be able to have all software automatically installed into
the working directory without needing admin rights on the machine.

Add the `lcdb` channel (one-time setup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This enables the `lcdb` conda channel so that additional dependencies not
included in bioconda can be installed::

    conda config --add channels lcdb

Create a new conda environment (one-time setup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This creates a top-level environment with Snakemake and other requirements. It
should be activated any time you'll be working with `lcdb-wf`. Here we're using
the name "lcdb-wf", but you can use anything::

    conda create -n lcdb-wf -y python=3 --file requirements.txt

Then activate the environment::

    source activate lcdb-wf

When you're done you can deactivate, though you might want to hold off on this
for now::

    source deactivate

Download example data (one-time setup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will download the example data to the directory ``data/``::

    python ci/get-data.py


Run the RNA-seq workflow with example data
------------------------------------------

With the `lcdb-wf` environment activated, change to the ``workflows/rnaseq``
directory, and run::

    snakemake --use-conda

Specify more than one core with ``-j``, e.g., ``-j 8`` to use more cores.

This workflow includes the ``workflows/references/Snakefile`` workflow. It
downloads genome sequence and reference files and builds indexes as necessary
(HISAT2 genome index, salmon transcriptome index, bowtie2 index for rRNA, GTF
file of gene annotations) and then carries on with the RNA-seq workflow.

These references are configured in the ``config/config.yaml`` file, which is
heavily commented to guide you on how to modify it for your own needs.

The RNA-seq workflow includes the standard mapping, counting, and differential
expression stages, as well as many quality-control steps. See :ref:`rnaseq` for more details.

The DAG of jobs looks like this:

.. image:: rnaseq.png

After the workflow runs, here are some useful Points of interest in the output:

    - ``data/rnaseq_samples/*``: sample-specific output. For example,
      individual BAMs can be found here
    - ``data/aggregation/multiqc.html``:  MultiQC report.
    - ``downstream/rnaseq.html``: Differential expression results.

Run the ChIP-seq workflow with example data
-------------------------------------------

With the `lcdb-wf` environment activated, change to the ``workflows/chipseq``
directory, and run::

    snakemake --use-conda

Specify more than one core with ``-j``, e.g., ``-j 8`` to use more cores.

Like the RNA-seq workflow, the ChIP-seq workflow includes the
``workflows/references/Snakemake`` workflow, so that genome fastas are
downloaded and indexes built as necessary.

The ChIP-seq workflow includes:

    - trimming reads with cutadapt
    - mapping reads with Bowtie2
    - FastQC on raw, trimmed, and aligned reads
    - Remove multimappers (samtools) and duplicates (Picard MarkDuplicates)
    - fastq_screen on multiple configured genomes to look for evidence of
      cross-contamination
    - QC aggregation using MultiQC, along with a custom table for library sizes
    - merging and re-deduplicating to correctly handle technical replicates
    - bigWigs created from unique, no-dups BAM files
    - deepTools plotFingerprint run on grouped IP and input for QC and
      evaluation of enrichment
    - peak-calling using macs2 and/or spp
    - conversion of BED files into bigBed (or bigNarrowPeak where possible)
    - track hub of bigWigs and bigBeds to visualize peak-calling in UCSC Genome Browser

.. image:: chipseq.png

Points of interest:

    - ``data/chipseq_samples/*``: sample-specific output. Individual BAM files
      for a sample can be found here.
    - ``data/chipseq_merged/*``: technical replicates merged and re-deduped, or
      if only one tech rep, symlinked to the BAM in the samples directory
    - ``data/chipseq_peaks/*``: peak-caller output, including BED files of
      called peaks and bedGraph files of signal as output by each algorithm
    - ``data/chipseq_aggregation/multiqc.html``: MultiQC report

Run the references workflow with example data
---------------------------------------------

This is optional; parts of this workflow were actually run automatically as
needed for the RNA-seq and ChIP-seq workflows. With the `lcdb-wf` environment
activated, change to the ``workflows/references`` directory and run::

    snakemake --use-conda

Adjust the ``-j`` argument to match the number of CPUs to run jobs in parallel
and speed up the workflow.

.. image:: references.png


"External" workflow
-------------------
This workflow is a working example downloads some data from modENCODE in an
older fly genome assembly (dm3), fixes the formatting so they can be lifted
over, and lifts over the files to the newer dm6 assembly.

It can be used as a template for integrative downstream work, as a place to
keep track and automate the download and preparation of external published
data. It can then be incorporated into the ``figures`` workflow (described
below) to integrate the analysis with other output.

.. image:: external.png


"Figures" workflow
------------------

This workflow is a working example of how you would tie together output from
RNA-seq, ChIP-seq, and "external" workflows to automate figure creation and
formally link together the dependencies of each figure all the way back to the
original fastq files. If any changes are made upstream, they will trigger
downstream rules to re-run as needed.

For example, if a new GTF annotation file comes out, you would change the URL
in the config and re-run the RNA-seq workflow. Only those jobs that depended on
some way on the GTF file will be re-run, up to and including any figures that
ultimately depended on that annotation file (so trimming and alignment will
not run again, but featureCounts, DESeq2, and any figures that use DESeq2
output will be re-run).

To provide a sufficiently complex example that can be used in real-world
applications, this workflow currently:

    - counts the number of peaks in all configured peak-calling runs and stores
      the output in a TSC report (work is performed in the script
      ``scripts/peak_count.py``)
    - identifies peaks at promoters of genes, reports a summary of how many
      peaks from each run were found in a promoter, and creates BED files of
      these subsets (``scripts/peaks_at_promoters.py``)
    - builds the DAG image for the ChIP-seq and RNA-seq workflows
    - symlinks over the ChIP-seq peaks and RNA-seq differential expression results
    - extracts the text from the README.txt files created by the scripts
      (usually from their docstrings), and compiles them into a summary report
    - the ``figures`` directory can then be zipped up and distributed to collaborators


Colocalization workflow
-----------------------
The output of this workflow is a set of heatmaps showing metrics of
colocalization between pairs of regions. This can be used to answer questions
like "what else does my protein of interest bind with?".

Several colocalization methods are run, and they all give slightly different
results.

- bedtools fisher
- bedtools jaccard
- GAT (Genome Association Test) log2 fold change
- GAT (Genome Association Test) nucleotide overlap
- IntervalStats

Next steps
----------
- Edit ``config/sampletable.tsv`` and ``config/config.yml`` to reflect your
  experimental design and desired references
- Add your original FASTQ files to ``data/rnaseq_samples/`` directories, likely
  using symlinks
- Edit the Snakefile for the workflow you're running to disable any steps you
  don't need
- After running the workflow, you can customize ``downstream/rnaseq.Rmd`` for
  your particular analysis.
