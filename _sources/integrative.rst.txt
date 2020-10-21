.. _integrative:

Integrative workflows
=====================

Here we look at integrative workflows which can be used to combine multiple
standard or non-standard workflows.

.. _colocalization:

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

.. _external:

"External" workflow
-------------------
Often we want to compare new data with existing published data. We have found
that in practice, having a separate workflow to handle downloading and
reformatting and various conversion tasks helps with organization.

The test workflow is a working example that:

- downloads ChIP-seq data from modENCODE in an older fly genome organism (dm3)
- downloads the chainfile for liftover
- fixes the formatting of the downloaded files so they can be lifted over
- lifts over the files to the newer dm6 assembly.

The file is intended to be heavily edited for the particular experiment; it is
here mostly as a placeholder and to be used as a template for integrative
downstream work.  It can then be incorporated into the ``figures`` workflow (see
:ref:`figures`) to integrate the analysis with other output.

.. image:: external.png

.. _figures:

"Figures" workflow
------------------

This workflow is a working example of how you would tie together output from
RNA-seq, ChIP-seq, and "external" workflows to automate figure creation and
formally link together the dependencies of each figure all the way back to the
original fastq files. If any changes are made upstream, they will trigger
downstream rules to re-run as needed.

For example, if a new GTF annotation file comes out, you would change the URL
in the config, and then re-run the figures workflow. All RNA-seq jobs that
depended on some way on the GTF file will be re-run. This will include the
feature counts and downstream RNA-seq analysis, but will *not* re-run any of
the jobs like trimming or mapping that do not depend on the GTF. In addition,
any figures that depended in some way on that GTF file will also be re-run.



To provide a sufficiently complex example that can be used in real-world
applications, this workflow currently:

    - counts the number of peaks in all configured peak-calling runs and stores
      the output in a TSV report (work is performed in the script
      ``scripts/peak_count.py``)
    - identifies peaks at promoters of genes, reports a summary of how many
      peaks from each run were found in a promoter, and creates BED files of
      these subsets (``scripts/peaks_at_promoters.py``)
    - builds the DAG image for the ChIP-seq and RNA-seq workflows
    - symlinks over the ChIP-seq peaks and RNA-seq differential expression results
    - extracts the text from the README.txt files created by the scripts
      (usually from their docstrings), and compiles them into a summary report
    - the ``figures`` directory can then be zipped up and distributed to collaborators

