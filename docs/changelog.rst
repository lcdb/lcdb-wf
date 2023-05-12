Changelog
=========

v1.10.1
-------
This is a bugfix and minor patch release.

- Bugfix: the references workflow was missing the ``resources:`` directives;
  they have now been added.

- The new ``utils.autobump`` function can be used to easily specify default and
  incremented resources, and the ``utils.gb`` and ``utils.hours`` make it
  a little easier to specify when autobump is not required.

  In the following example, memory will be set to 8 * 1024 MB and will
  increment by that much each retry. The runtime will be set to 2 * 60 minutes,
  and will increment by 10 * 60 minutes each retry. The disk will be set to 100
  * 1024 MB, and will not increase each retry.

  .. code-block:: python

      resources:
          mem_mb=autobump(gb=8),
          runtime=autobump(hours=2, increment_hours=10),
          disk_mb=gb(100)

- WRAPPER_SLURM no longer has the ``--latency-wait=300``,
  ``--max-jobs-per-second=1``, and ``--max-status-checks-per-second=0.01``
  which would override any profile settings.

- In RNA-seq and ChIP-seq, the cutadpt rule now defaults to using
  ``--nextseq-trim 20`` instead of ``-q 20``, to better handle the majority of
  sequencing data we have recently been working with (NovaSeq). See `this
  section of the cutadapt docs
  <https://cutadapt.readthedocs.io/en/stable/guide.html#nextseq-trim>`_ for
  details.

v1.10
-----
The major change here is refactoring the Snakefiles to use the ``resources:``
directive in each rule, and removing the ``--clusterconfig`` mechanism which
has long been deprecated.

For running on a cluster, this requires a `profile
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_.
E.g., on `NIH's Biowulf <https://hpc.ni.gov>`_, use the `NIH-HPC
snakemake_profile <https://github.com/NIH-HPC/snakemake_profile>`_.

General
~~~~~~~
- No longer using clusterconfig, instead using resources to configure cluster resources
- Migrated to a unified testing script that simplifies local and CI testing
- If sampletable is from SRA, raise an error if a Layout column can't be found
  (to prevent incorrect interpretation of samples as single-end)
- Ensure bam indexes are made for the markdups bams, even if bigwigs are not created
- Remove libsizes table, which was largely redundant with fastqc results

RNA-seq
~~~~~~~
- Fix R tests
- All ``lcdbwf`` R functions use the ``:::`` namespace lookup syntax
- Fix library loads in rnaseq.Rmd to ensure they come before parallelization configuration
- New function ``lcdbwf:::lfc_scatter`` for comparing multiple DESeq2 contrasts
- Updates and fixes to ``gene-patterns.Rmd``


v1.9
----

This version has substantial changes in the ``rnaseq.Rmd`` file to streamline
its use in a production environment. This involves moving most of the code
complexity into the ``lcdbwf`` R package and using a new config file as much as
possible. See details below.

General
~~~~~~~
- environments have been updated with recent versions of all tools
- WRAPPER_SLURM arguments updated with arguments better suited for cluster submission
- PhiX reference configs have been removed
- compatibility with Python 3.10
- fastq-dump rules have been converted to scripts. This is because sra-tools in
  versions earlier than 3.0 have issue with SSL certs, however sra-tools=3
  cannot be installed alongside recent versions of salmon (due to conflicting
  pinnings with the ``icu`` package). Therefore, fastq-dump is now run as
  a script in its own conda environment.
- new idxstats rule for chipseq and rnaseq

RNA-seq
~~~~~~~

**This version has major changes to** ``rnaseq.Rmd``. Briefly:

1. This file has been overhauled to be driven by a config file. This
   dramatically reduces the need to scroll through the RMarkdown file and make
   all the customizations for a particular experiment. Now, editing the config
   file sets up most of the project-specific components. Note that contrasts
   still need to be customized in the Rmd file.
2. The narrative and explanatory text has been moved to ``text.yaml`` and is
   included at render time. This reduces the need to scroll through lots of
   boilerplate text in the RMarkdown while still retaining the ability to
   easily edit it.
3. Most of the complexity has been offloaded to the ``lcdbwf`` R package.
4. Caches are much improved. See the :ref:`downstream-detailed` section for
   more information.
5. Functional enrichment is moved into a separate RMarkdown file.

Downstream RNA-seq config
,,,,,,,,,,,,,,,,,,,,,,,,,

The file, `workflows/rnaseq/downstream/config.yaml` is heavily commented to
describe the various settings. The sections of the config are designed such
that they can be used as additional chunk options to chunks in which they are
used. This additional chunk option is used by RMarkdown to compute the hash of
the chunk. The result is that making a change in the config file is sufficient
to invalidate the cache of any chunks that specify that section as a chunk
option.

Complexity moved to ``lib/lcdbwf/R``
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

Another major change is that most of the complexity in the ``rnaseq.Rmd`` file
has been factored out into the ``lcdbwf`` R package that is stored inn
``lib/lcdbwf``. While this means that all code is no longer included in the
final rendered HTML file, it does make the Rmd much more streamlined to work
with. It also has the side effect of making it easier to write unit tests on
separate functions.

Many helper functions have been added to the ``lcdbwf`` R package, including
ones to streamline the creation of dds and results objects, composing and saving
them, and generating many of the outputs.

Improved caching of results chunks
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

A somewhat major change is a new strategy for allowing ``results()`` calls to be
split across multiple, independently-cached chunks that are then properly merged
together into a single ``res.list`` object while handling dependencies and
parallelization (thanks to `@njohnso6 <https://github.com/njohnso6>`_). This
dramatically speeds up the process of incrementally adding contrasts to complex
experimental designs.

Other changes
,,,,,,,,,,,,,

In addition to these major changes, there are also many other improvements
to ``rnaseq.Rmd``:

    - AnnotationHub databases are only retrieved from cache when they are
      needed. This dramatically speeds up rendering of the HTML, since before
      the OrgDb would always load no matter what.
    - Toggle Kallisto or Salmon quantification with a simple true/false; this
      automatically sums to gene level using automatically retrieved TxDb. This
      also now supports creating dds objects from featureCounts, Salmon, or
      Kallisto in such a way that they can be easily compared with each other.
    - ``lcdbwf::compose_results()`` to combine res_list and dds_list objects
      together by inspecting the global namespace for specially-named objects
    - Helper functions for retrieving global config and data structures (e.g.,
      ``lcdbwf::get_config()``, ``lcdbwf::get_dds()``)
    - Helper function ``lcdbwf::match_from_dots`` for working with `...`
      arguments and splitting them up to only go to the functions they are
      intended for
    - Much faster to attach info (e.g., adding SYMBOL to all results) since the
      AnnotationDbi calls are only done once instead of for each results
      object.
    - Refactored functional enrichment to be much more generalized, currently
      using Gene Ontology and MSigDB. MSigDb, via the ``msigdbr`` package, is
      available for multiple species and so this incorporates Reactome and
      KEGG. But the generalized method can be applied to any arbitrary gene
      sets, allowing for much more customization.
    - Fixes to clusterProfiler::emapplot calls in particular corner cases
    - Functional enrichment is now a completely separate file, using the
      ``combined.Rds`` file as an intermediate between ``rnaseq.Rmd`` and
      ``functional_enrichment.Rmd``.
    - All-in-one enrichment function that runs either overrepresentation or
      GSEA. Makes it much easier to do *ad hoc* tests.
    - Helper function ``lcdbwf::enrich_list_lapply()`` to apply arbitrary
      functions to the highly-nested `enrich_list` data structure
    - Helper function ``lcdbwf::collect_objects`` to help compile discovered
      results objects
    - ``lcdbwf::get_sig()`` has more options for what to return
    - Plotting wrappers for clusterProfiler plot functions, allowing plots to be
      configured via the config file.
    - New dds diagnostics and results diagnostics functions and sections of the
      Rmd, useful for troubleshooting
    - Refactored the results tabs: MA plots come first; ensure 10 genes are always plotted in MA
      plots, added volcano plots with labeled genes, removed top 3 and bottom
      3 gene plots
    - PCA plots using plotly no longer need "unrolled" for-loops; multiple PCA
      coloring and clustered heatmap row side colors are now configured in the
      YAML config file
    - Moved size factor plots and gene version removal to lcdbwf package
    - Use datatable to show initial sampletable for cleaner output
    - Make original dds_initial object the same way as later dds objects and
      always using a design of ``~1`` to be used in PCA and heatmaps
    - "Differential expression" header moved so that code is no longer hidden
      under the size factors plot
    - Option for filling in NA in symbol with Ensembl IDs
    - collapseReplicates2 uses ``collapse_by`` rather than ``combine.by``
    - Updated the code style throughout to use the tidyverse/google style guide
    - RNA-seq differential expression output is additionally included in an
      Excel file with one sheet per contrast.

Tests
~~~~~

- ``lcdbwf`` R package now has its own tests via ``devtools`` and ``testthat``
- recent versions of Snakemake are broken when ``--until`` is used in certain
  circumstances; a ChIP-seq test has been disabled temporarily.
- after a successful test, the environment is written as an artifact on circleci

References
~~~~~~~~~~

- Fixed a longstanding issue with *S. cerevisiae*, now the GFF file is properly converted to GTF.

v1.8
----

General
~~~~~~~

- Complete shift to using pinned ``env.yaml`` files to specify conda
  environments, and using ``mamba`` for building environments (consistent with
  recent versions of Snakemake). This is now reflected in documentation and
  the updated-and-improved ``deploy.py``.

- Reorganization/cleanup of the ``include`` directory

- Added conda troubleshooting notes to the documentation (see
  :ref:`conda-troubleshooting`).

- The ``lib.helpers.preflight`` function no requires the first column of the
  sampletable to be named `samplename` when checking configs.

- Improvements to the deployment script ``deploy.py``:

    - now requires Python >3.6
    - proper logs (so you can easily see how long it takes to build an env)
    - supports downloading and running the script directly, which will clone
      a temporary copy and deploy from there
    - using Control-C to stop the deployment will also stop mamba/conda
    - colored output
    - mamba is used by default, but ``--conda-frontend`` will use conda instead

- fastq-dump log is sent to file rather than printed to stdout

- Threads: cutadapt single-end now uses specified threads (it was using
  1 thread by default); use 6 threads for fastqc

- Added new preflight checks for RNA-seq and ChIP-seq specific configs.

- Added a ``run_complex_test.sh`` driver script for testing the workflows on
  full-scale publicly available data 

RNA-seq
~~~~~~~

- **Configuration change:** The ``stranded:`` field is now required for RNA-seq.
  This is used to choose the correct parameters for various rules, and avoids
  one of the main reasons to edit the Snakefile. See :ref:`cfg-stranded` for
  more details on its use.

- added ``stranded:`` field to all configs used in testing

- The ``strand_check`` rule now runs MultiQC for a convenient way of evaluating
  strandedness of a library.

- Kallisto is now supported in both the RNA-seq Snakefile, references
  Snakefile, included reference configs, and downstream ``rnaseq.Rmd``


References
~~~~~~~~~~

- When checking URLs in reference configs, don't use ``curl`` to check
  ``file://`` URIs.

- There is a new feature for reference configs that allows chaining
  post-processing functions together, see :ref:`advanced-postprocessing`. This
  means that it is possible, for example, to add ERCC spike-ins (which need
  post-processing) onto references that themselves need post-processing.

- ``lib/postprocess/ercc.py`` has new helper functions for adding ERCC
  spike-ins to fasta files and GTF files.

- added ``'kallisto'`` to included reference configs

ChIP-seq
~~~~~~~~

- symlinks rule is now local
- added collectinsertsizes pattern to support PE ChIP-seq experiments
- merging bigwigs log no longer goes to stdout


v1.7
----

Setup
~~~~~

Use mamba for installation of environments, consistent with Snakemake recommendations

Testing
~~~~~~~

- We now recommend using `mamba <https://github.com/mamba-org/mamba>`_ to
  create conda environments. This is dramatically faster and solves some
  dependency issues. Our automated tests now use this.

- We have moved from requirements.txt files to env.yaml files. We also now
  encourage the use of the strictly-pinned environments for a more stable
  experience to hopefully avoid transient issues in the packaging ecosystem.

- ``tbb=2020.2`` as a dependency to fix a recent packaging issue with conda-forge.

- many documentation improvements

- symlinks rule is only set to localrule when it exists (it does not exist when
  running an analysis exclusively from SRA)

References
~~~~~~~~~~

- updated URLs for those that have changes (e.g., Sanger -> EBI; using https
  instead of ftp for UCSC-hosted genomes)

- new ``gff2gtf`` post-process tool for when an annotation is only available as
  GFF. *S. pombe* needs this, for example, and the
  `Schizosaccharomyces_pombe.yaml`` reference config has been updated
  accordingly.


- The references workflow no longer reads the config file in its directory.
  This fixes some subtle overwriting issues when providing config files or
  items from the command line during as is used during certain test runs. If
  running the references workflow alone, it must be called with
  ``--configfile``

RNA-seq
~~~~~~~

- featureCounts now uses BAM files with duplicates marked. Previously if you
  wanted to run featureCounts in a mode where it excluded duplicates you would
  need to reconfigure rules.

- improved comments in RNA-seq downstream RMarkdown files

Testing
~~~~~~~

- new test that checks all URLs identified in config files to ensure that the
  included reference files remain valid

- there is now a separate ``run_downstream_test`` script`

- simplified the CircleCI DAG to optimize testing resources

v1.6
----

References
~~~~~~~~~~
- overhaul the way transcriptome fastas are created. Instead of requiring
  separate download, they are now created out of the provided GTF and fasta
  files. The reference config section now uses keys ``genome:``,
  ``transcriptome:``, and ``annotation:`` rather than the ``fasta:`` and
  ``gtf:`` keys.
- **backwards-incompatible change:** reference config files have been updated
  to reflect the changes in the references workflow
- Update PhiX genome fasta to use NCBI rather than Illumina iGenomes

ChIP-seq workflow
~~~~~~~~~~~~~~~~~
- ChIP-seq workflow now properly supports paired-end reads
- A ChIP-seq workflow can now be run when the ``chipseq:`` and/or
  ``peak_calling:`` sections are omitted.
- added a missing bowtie2 config entry in ``clusterconfig.yaml`` which could
  result in out-of-memory errors when submitting to a cluster using that file


RNA-seq workflow
~~~~~~~~~~~~~~~~
- if colData is a tibble this no longer causes issues for importing counts
- dupRadar removed from RNA-seq workflow. We ended up never using it, and it
  depends on R which we've since removed from the main environment.
- new ``strand_test`` rule, which can be run explicitly with ``snakemake -j2
  strand_check``. This generates ``strandedness.tsv`` in the current directory,
  which is the summarize output of RSeQC's ``infer_experiment.py`` across all
  samples.
- implement STAR two-pass alignment. Default is still single-pass.
- Clean up hard-coded STAR indexing Log.out file
- Include ``ashr`` and ``ihw`` Bioconductor packages in the R requirements, for
  use with recent versions of DESeq2.


RNA-seq downstream
~~~~~~~~~~~~~~~~~~

- Functional enrichment and gene patterns are now separate child documents.
  This makes it easier to turn them on/off by only needing to adjust the chunk
  options of the child chunk
- Created a new documentation method for rnaseq.Rmd. Now there is a separate,
  dedicated documentation page with sections that exactly correspond to each
  named chunk in the Rmd, as well as a tool for ensuring that chunks and docs
  stay synchronized. See :ref:`rnaseqrmd` for the new docs.
- New ``counts.df`` and ``counts.plot`` functions to make it much easier to
  make custom dotplots of top counts by melting and joining the counts table
  with the metadata in colData.
- DEGpatterns cluster IDs are now added as additional columns in the output
  TSVs for each contrast
- Many functions in the rnaseq.Rmd now expect a list of :term:`dds` objects.
  See :ref:`dds_list` for more info on this.
- Created a new R package, ``lcdbwf`` stored in :file:`lib/lcdbwf`. This can be
  edited in place, and it is loaded from disk within ``rnaseq.Rmd``.
- Modified some output keys to support recent versions of Snakemake, for which
  ``count`` is a reserved keyword


General
~~~~~~~
- Conda environments are now split into R and non-R. See :ref:`conda-envs` for
  details. Updated ``deploy.py`` accordingly
- symlinks rules are now set to be localrules
- updated workflows to work on recent Snakemake versions
- split environments into non-R and R. This, along with a loose pinning of
  versions (``>=``), dramatically speeds up environment creation.
- updates to support latest Snakemake versions
- improvements to testing:
   - environment YAML files, rendered HTML, and docs are stored as artifacts on CircleCI
   - consolidations of some RNA-seq tests to reduce total time
   - additional comments in the test config yaml to help new users understand the system
- new "preflight check" function is run to hopefully catch errors before running workflows
- updates to support recent Picard versions
- added wildcard constraints to help Snakemake solve DAG


v1.5.3
------

General
~~~~~~~
- default 12-hr wall time in WRAPPER_SLURM
- update .gitignore (`#223 <https://github.com/lcdb/lcdb-wf/issues/223>`_)
- remove the FastQC status checks section from the MultiQC report (which shows
  up in recent MultiQC versions) (`#246 <https://github.com/lcdb/lcdb-wf/issues/246>`_

Bugs
~~~~

- add bed12 conversion for all species with default reference configs
- presence of an orig_filename_R2 in sampletable is sufficient to consider the
  experiment PE
- ensure DEGpattern output only contains unique genes
- bring back featurecounts in multiqc report
- "attach" chunk in rnaseq.Rmd was not properly set to depend on the "results" chunk

RNA-seq
~~~~~~~

- dds objects can now be created from a full featureCounts input file and
  a subsetted colData table, if subset.counts=TRUE
- improve the dependencies between rnaseq.Rmd chunks so that cache=TRUE behaves
  as expected: (`#232 <https://github.com/lcdb/lcdb-wf/issues/232>`_)
- add plots for rnaseq.Rmd size factors (`#222 <https://github.com/lcdb/lcdb-wf/issues/222>`_)
- run rseqc instead of CollectRnaSeqMetrics (the multiqc output is nicer for
  it, and it's pretty much doing the same thing) (`#218 <https://github.com/lcdb/lcdb-wf/issues/218>`_)
- when converting Ensembl to symbol, if there is no symbol then fall back to
  the Ensembl ID to avoid NA (`#246
  <https://github.com/lcdb/lcdb-wf/issues/246>`_)
- in rnaseq.Rmd, all caches will be invalidated if the sampletable or the
  featurecounts table have changed.

Tests
~~~~~
- using continuumio/miniconda3 container; finally got en_US.utf8 locale
  installed and working correctly in that container so that multiqc works.


v1.5.2
------

Bug fixes
~~~~~~~~~

- When some samples were substrings of other samples (e.g., `WT_1_1` and
  `WT_1_10`), DESeqDataSetFromCombinedFeatureCounts was assigning the wrong
  names. This has now been fixed in `helpers.Rmd`.

v1.5.1
------

Bug fixes
~~~~~~~~~

- DESeqDataSetFromCombinedFeatureCounts (added in v1.5) was incorrectly
  assigning labels to samples when the order of the sampletable did not match
  the order of the samples in the featureCounts table columns. This has been
  fixed.

General
~~~~~~~

- `deploy.py` deployment script now only pays attention to files checked in to
  version control and optionally can create a conda environment in the target
  directory.

- tests now work out of a newly-deployed instance to better reflect real-world
  usage


ChIP-seq and RNA-seq
~~~~~~~~~~~~~~~~~~~~
- reorder cutadapt commands to avoid a MultQC parsing bug in the cutadapt
  module (see https://github.com/ewels/MultiQC/issues/949)

RNA-seq
~~~~~~~
The majority of these changes affect ``rnaseq.Rmd``:

- modifications to MultiQC config to get back featureCounts output
- `plotMA.label` function (in ``helpers.Rmd``) now defaults to FDR < 0.1
  (instead of 0.01), and additionally supports labeling using different columns
  of the results object (e.g., "symbol").
- remove some now-redundant featureCounts code
- add a comment showing where to collapse replicates
- convert colData's first column to rownames
- implement lower limit for DEGpatterns clustering (default is 0, but can
  easily set to higher if you're getting issues)
- expose arbitrary additional function arguments to ``top.plots``. This allows
  different `intgroup` arguments to be passed to the `my.counts` function,
  enabling different ways of plotting the gene dotplots.


v1.5 (Sept 2019)
----------------

Major change: **it is no longer possible to mix single-end and paired-end
samples within the same run of the workflow.** See `#208
<https://github.com/lcdb/lcdb-wf/pull/208>`_ and the corresponding issue
description at `#175 <https://github.com/lcdb/lcdb-wf/issues/175>`_.

This version also has many improvements to the ``rnaseq.Rmd`` file for RNA-seq,
as described below.

RNA-seq
~~~~~~~

Many changes and improvements to ``rnaseq.Rmd``, including:

- Differential analysis summaries now include labeled MA plots (`#192 <https://github.com/lcdb/lcdb-wf/pull/192/files>`_)
- PCA plots now use plotly for improved insepction of samples (`#192 <https://github.com/lcdb/lcdb-wf/pull/192/files>`_
- don't use knitrBootstrap any more (`#192 <https://github.com/lcdb/lcdb-wf/pull/192/files>`_
- heatmaps use heatmaply package for better interaction (`#192 <https://github.com/lcdb/lcdb-wf/pull/192/files>`_
- allow ``sel.list`` to be used for UpSet plots and fix some typos `#205 <https://github.com/lcdb/lcdb-wf/pull/205>`_
- workaround for degPatterns for corner cases where there are few clusters because of the ``minc`` parameter (`#205 <https://github.com/lcdb/lcdb-wf/pull/205>`_)
- alpha and lfc.thresh are now pulled out into a separate chunk (`#206 <https://github.com/lcdb/lcdb-wf/pull/206>`_)
- Support AnnotationHub http proxy handling in new version of AnnotationHub (`#207 <https://github.com/lcdb/lcdb-wf/pull/207>`_).

As well as the following changes to other parts of the RNA-seq workflow, such as:

- better bigWig file nomenclature (`#194 <https://github.com/lcdb/lcdb-wf/pull/194/files>`_), uses "pos" and "neg".
- featureCounts only runs once on all BAMs rather than individual samples (`#195 <https://github.com/lcdb/lcdb-wf/pull/195>`_)
- support `rseqc infer_experiment`, which replaces running featureCounts in multiple stranded modes (`#199 <https://github.com/lcdb/lcdb-wf/pull/199>`_, `#203 <https://github.com/lcdb/lcdb-wf/pull/203>`_)
- use ``--validateMappings`` for salmon (`#203 <https://github.com/lcdb/lcdb-wf/pull/203>`_)

References
~~~~~~~~~~
- fix typo in *S. pombe* name

All workflows
~~~~~~~~~~~~~

- Documentation now recommends creating an environment for each directory using the `-p` argument (`#195 <https://github.com/lcdb/lcdb-wf/pull/195>`_)


v1.4.2 (Jul 2019)
-----------------

Bugfixes
~~~~~~~~

- Don't require ChIP-seq configs to have at least one block for each supported
  peak-caller

v1.4.1 (Jul 2019)
-----------------

RNA-seq
~~~~~~~

- KEGG results were not being added to the ``all.enrich`` list in ``rnaseq.Rmd``
- symlinking bigWigs is now a local rule
- default cutadapt options have changed to reflect current recommendations from
  the author, and the cutadapt rule is now explicity using arguments rather
  than requiring a separate ``adapters.fa`` file.
- featureCounts now auto-detects whether it should be run with the ``-p``
  argument in paired-end mode (previously it was up to the user to make sure
  this was added). The rule does have an override if this behavior is not wanted.

References
~~~~~~~~~~

- The reference config for *Drosophila* is now fixed. Previously it depended on
  `chrom_convert`. That script was a fly-specific script in lcdblib, but
  lcdblib is no longer a dependency since v1.3. This fix uses the
  `convert_fastq_chroms` and `convert_gtf_chroms` used in reference configs for
  other species.

v1.4 (May 2019)
---------------
RNA-seq
~~~~~~~
Much-improved ``rnaseq.Rmd``:

- tabbed PCA plot
- improved DEGpatterns chunk
- dramatically improved functional enrichment section, with tabbed clusterprofiler plots and exported data in two flavors (combined and split)
- improved upset plots, with exported files showing sets of genes
- improved comments to highlight where to make changes
- add new helper functions to ``helpers.R``:
   - ``fromList.with.names``, for getting UpSet plot output
   - ``rownames.first.col``, to make tidier dataframes
   - ``nested.lapply``, for convenient 2-level nested list apply
   - clusterprofiler helper functions


v1.3 (May 2019)
---------------
Bugfixes
~~~~~~~~
- Fix broken paired-end support for RNA-seq. Previously, when using data from
  elsewhere on disk and using the symlink rules, R2 would be symlinked to the
  same file as R1.
- Support for Snakemake 5.4.0 which changes behavior of the ``expand()``
  function.

Infrastructure
~~~~~~~~~~~~~~
- new deploy script to copy over only the files necessary for an analysis,
  avoiding the clutter of testing infrastructure.
- lcdblib, an external package, is no longer a dependency. In the interest of
  better transparency and to make the code here easier to follow, the relevant
  code from lcdblib was copied over to the ``lib`` directory in this
  repository.

ChIP-seq and RNA-seq
~~~~~~~~~~~~~~~~~~~~

- Bowtie2, HISAT2, and rRNA rules no longer use wrappers. This makes it easier
  to track down what parameters are being used in each rule.
- RSeQC is now available in Python 3 so wrappers have been removed.
- NextGenMap support removed

v1.2 (Mar 2019)
---------------

RNA-seq
~~~~~~~
- First-class paired-end support, including mixing PE and SE samples in the
  same sampletable

- Support for STAR aligner

References
~~~~~~~~~~
- FASTA files are always symlinked into the directories of indexes that were
  created from it

- Reference configs:

   - updated existing
   - added more species
   - new post-process for fasta or gtf: you can now use
     NICHD-BSPC/chrom-name-mappings to convert chromosome names between UCSC
     and Ensembl (see reference configs for examples of use)

ChIP-seq and RNA-seq
~~~~~~~~~~~~~~~~~~~~
- Updates to dependencies and MultiQC config

Infrastructure
~~~~~~~~~~~~~~

- Updated requirements in ``requirements.txt`` and in wrappers

- Changed all ``pd.read_table()`` to ``pd.read_csv(sep="\t")`` to prevent warnings

- Changed all ``yaml.load()`` to ``yaml.load(Loader=yaml.FullLoader)`` to
  prevent warnings

- Using DeprecationWarning rather than UserWarning in the deprecation handler
  so there's less spam in the logs

- Improved tests:

  - using data from pybedtools repo because modENCODE seems to be down
  - append rather than prepend base conda to PATH on circleci
  - separate isolated tests for STAR, ngm, and SRA
  - updated conda

- Docs additions:

  - TMPDIR handling
  - clusterconfig
  - WRAPPER_SLURM
  - docs for developers
  - symlinking fastqs
  - using SRA sampletables
  - paired-end data

Colocalization
~~~~~~~~~~~~~~
- From colocalization, removed the GAT "fractions" heatmap due to unresolved
  pandas index errors

v1.1 (Aug 2018)
---------------

Infrastructure
~~~~~~~~~~~~~~

- The default settings in Snakefiles are for real-world use, rather than for
  testing. This reduces the amount of editing necessary before running actual
  data. See :ref:`test-settings` for the extra step to take when testing
  locally.

- new ``run_test.sh`` script in each workflow directory to automatically run
  the preprocessor when running test data

- added extensive comments to Snakefiles with ``NOTE:`` string to make it
  obvious where and how to make changes.

- Documentation overhaul to bring everything up to v1.1. This includes Sphinx
  autodocs on the ``lib`` module.

- pytest test suite is run on the ``lib`` module

References
~~~~~~~~~~

- new `metadata` section in references config, which can be used to store
  additional information like mappable bases and genome size.

- References can now be included from other YAML files into the main config
  file. This dramatically simplifies individual configfiles, and allows
  multiple workflows to use identical references without having to do
  error-prone and hard-to-maintain copy/pastes between workflow configs. See
  :ref:`references-config` for details.

- New GTF conversion, ``mappings``. This is intended to replace the
  ``annotation_hub`` conversion, which was problematic because 1) a particular
  annotation hub accession is not guaranteed to be found in new versions of
  AnnotationHub, resulting in lack of reproducibility, and 2) it was difficult
  to synchronize the results with a particular GTF annotation. The
  ``annotation_hub`` conversion is still supported, but if it's used then
  a DeprecationWarning will be emitted, recommending ``mappings`` instead.


Both RNA-seq and ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~

- `fastq_screen` is now configured via ``config.yaml``. This reduces the need
  to edit the Snakefile and coordinate between the config and the fastq_screen
  rule. Now everything is done within the config file.

- `fastq_screen` wrapper now handles additional output files created when using
  the ``--tag`` and ``--filter`` arguments to ``fastq_screen``.

- In the config file, ``assembly`` has been changed to the more-descriptive
  ``organism``. The change is backwards compatible, but a DeprecationWarning is
  raised if ``assembly:`` is still used, and changed to ``organism`` (though
  only in memory, not on disk).

- Patterns no longer use ``{sample_dir}``, ``{agg_dir}``, etc placeholders that
  need to be configured in the config YAML. Instead, these directories are
  hard-coded directly into the patterns. This simplifies the config files,
  simplifies the patterns, and removes one layer of disconnect between the
  filenames and how they are determined.

- removed 4C workflow since it used 4c-ker

ChIP-seq
~~~~~~~~
- macs2 and sicer can accept mappable genome size overrides

RNA-seq
~~~~~~~

- RNA-seq downstream:

    - ``downstream/help_docs.Rmd`` can be included for first-time users to
      describe the sections of the RNA-seq analysis

    - ``rnaseq.Rmd`` now uses the same ``NOTE:`` syntax as the Snakefiles for
      indicating where/what to change

    - Easy swapping of which strand to use from the three featureCounts runs
      performed by the workflow

    - Be explicit about using DESeq2::lfcShrink as is now the default in recent
      DESeq2 versions

    - improved the mechanism for keeping together results objects, dds objects, and
      labels (list of lists, rather than individual list object; refactored
      functions to use this new structure

v1.0.1 (Jun 2018)
-----------------
Bugfixes, last release before references changes.

Infrastructure
~~~~~~~~~~~~~~

- Transition to CircleCI for testing
- Use production settings by default; see :ref:`test-settings` for
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
- featureCounts is now run in all three strandedness modes, and results
  incorporated into MultiQC as separate modules.
- RNA-seq now symlinks "pos" and "neg" bigWigs, which describe how reads map to
  the *reference*, to "sense" and "antisense" bigWigs, which describe the
  *originating RNA*. This makes it easy to swap strands depending on protocol.
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

v1.0 (May 2018)
---------------
First official full release.
