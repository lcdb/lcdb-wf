Changelog
=========

v1.0
----
First full release.

v1.1
----
- ``include/reference_configs`` contains example config files for various
  species that can be pasted into individual workflows' ``references:``
  section.

- `fastq_screen` is now configured via ``config.yaml``. This reduces the need
  to edit the Snakefile and coordinate between the config and the fastq_screen
  rule. Now everything is done within the config file.

- In the config file, ``assembly`` is now ``organism``. The change is backwards
  compatible, but a DeprecationWarning is raised if ``assembly:`` is still
  used, and it is then automatically copied over to ``organism``

- using CircleCI for continuous integration testing

- the default settings in Snakefiles are for real-world use, rather than for
  testing. This reduces the amount of editing necessary before running actual
  data. See :ref:`test-settings` for the extra step to take when testing
  locally.

- All Snakefiles now have the string ``# NOTE`` to indicate settings that may
  need to be changed depending on the experiment.

- Documentation overhaul to bring everything up to v1.1, including Sphinx
  autodocs on the ``lib`` module.

- patterns no longer use ``{sample_dir}``, ``{agg_dir}``, etc placeholders that
  need to be configured in the config YAML. Instead, these directories are
  hard-coded directly into the patterns. This simplifies the config files,
  simplifies the patterns, and removes one layer of disconnect between the
  filenames and how they are determined.

- additional section in configs for averaging bigWigs in a flexible manner

- removed 4C workflow since it used 4c-ker

- ChIP-seq and RNA-seq can now run directly on sample tables from SRA's run
  browser (fastqs will be downloaded with fastq-dump, etc).

- pinning to deepTools >3.0 to use new command line args

- ChIP-seq now supports SICER for peak-calling.

- ChIP-seq multibigwigsummary and heatmaps for QC and clustering of samples

- ChIP-seq track hub with peaks and signal

- Peak caller wrappers now handle their own conversion to bigBed

- macs2 and sicer can accept mappable genome size overrides
  
- Improvements to MultiQC config so that multiple module runs (e.g., FastQC on
  raw, trimmed, aligned) show up separately in the stats table.

- Added colocalization workflow, external workflow, and figures workflow.

- RNA-seq now runs all three strandedness for featureCounts, and reports them
  all in MultiQC so you can check protocol expectations.

- RNA-seq now runs preseq for QC

- RNA-seq now symlinks "pos" and "neg" bigWigs, which describe how reads map to
  the *reference* to "sense" and "antisense" bigWigs, which describe the
  originating RNA. This makes it easy to swap strands depending on protocol.

- RNA-seq track hub with signal

- RNA-seq downstream:

    - ``downstream/help_docs.Rmd`` can be included for first-time users to
      describe the sections of the RNA-seq analysis

    - factored out ``downstream/helpers.Rmd`` into a separate file that is
      included into the main Rmd.

    - AnnotationHub is cached in ``include`` dir of repo to avoid clobbering
      any locally-installed versions.

    - UpSet plot example to show how to compare results from multiple contrasts

    - Easy swapping of which strand to use from the three featureCounts runs
      performed by the workflow

    - Be explicit about using DESeq2::lfcShrink as is now the default in recent
      DESeq2 versions

    - DEGpatterns plots
