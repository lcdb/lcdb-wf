Changelog
=========

v1.0
----
First full release.

v1.0.1
------
Bugfixes, last release before references changes.

Infrastructure
~~~~~~~~~~~~~~

- Transition to CircleCI for testing
- Use production settings by default; see :ref:`note-on-test-settings` for
  more.
- lots o' docs
- new ``include/references_configs`` to help organize references. These are
  currently not used by the workflows directly.
- bugfix: use additional options when uncompressing downloaded reference files
  (``--no-same-owner`` for ``tar``, ``-f`` for ``gunzip``)
- additional dependencies in the top-level environment to support the
  additional features in rnaseq.Rmd and track hubs.
- colocalization workflow, external workflow, figures workflow to demonstrate
  vertical integration

RNA-seq
~~~~~~~
- remove kallisto indexing, use salmon
- improvements to how chipseq sampletables are parsed (with more informative
  error messages)
- run preseq for RNA-seq library complexity QC
- support for merging bigwigs
- rule for symlinking stranded bigwigs
- featureCounts is now run in all three strandedness modes, and results
  incorporated into MultiQC as separate modules.
- new ``downstream/helpers.Rmd`` which factors out a lot of the work previously
  done in ``rnaseq.Rmd`` into separate functions.
- track hub building respects new sense/antisense bigwig symlinks

``downstream/rnaseq.Rmd``
~~~~~~~~~~~~~~~~~~~~~~~~~
- AnnotationHub uses cache dir that will not clobber default home directory cache
- use varianceStabilizingTransform instead of rlog
- print a size factors table
- use multiple cores for computationally expensive DESeq2 operations
- using separate lists for results, dds objects, and nice labels for automated
  plots for each contrast
- UpSet plots for comparing gene lists across contrasts
- DEGpattern plots for showing clusters of expression patterns (from the
  DEGreport package)
- attach normalized counts per sample and per factor (parsed from the model
  used for the contrast) as well as TPM estimates to the results tables
- trim the labels in GO enrichment plots when too long


ChIP-seq
~~~~~~~~
- sicer for chipseq domain calling
- pin snakemake <4.5.0 so that subworkflows behave correctly
- chipseq peak-calling rules (and therefore wrappers) now expect a chromsizes
  file as input
- bigbed files for narrowPeak and broadPeak files are created correctly
  depending on their format
- run multiBigWigSummary and plotCorrelation from deepTools for ChIP-seq QC
- ChIP-seq track hub generation script

Both RNA-seq and ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~
- update deeptools calls to reflect >v3.0 syntax
- support for SRA run tables so it's trivial to re-run experiments
  in SRA
- multiple FastQC runs are shown separately in MultiQC output
