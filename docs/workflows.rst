.. _workflows:

Overview of workflows
=====================

.. note::

   These workflows **are intended to be edited and customized by the user**.

   See :ref:`getting-started` for recommendations on setting up these workflows in
   your project directory.

Each workflow lives in its own directory:

::

    ├── references/
    │   ├── Snakefile
    │   └── ...
    ├── rnaseq/
    │   ├── Snakefile
    │   └── ...
    ├── chipseq/
    │   ├── Snakefile
    │   └── ...
    ├── colocalization/
    │   ├── Snakefile
    │   └── ...
    ├── external/
    │   ├── Snakefile
    │   └── ...
    └── figures/
        ├── Snakefile
        └── ...


There are two general classes of workflows, **primary analysis** and the
**integrative analysis**. 

Each workflow is driven by a ``Snakefile`` and is configured by plain text
`YAML <https://en.wikipedia.org/wiki/YAML>`_ and `TSV
<https://en.wikipedia.org/wiki/Tab-separated_values>`_ format files (see
:ref:`config` for much more on this).

Primary analysis workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~
The primary analysis workflows are generally used for transforming raw data
(fastq files) into usable results. For RNA-seq, that's differentially-expressed
genes (along with comprehensive QC and analysis). For ChIP-seq, that's called
peaks or differentially bound chromatin regions.

The primary analysis workflows are:

.. toctree::
   :maxdepth: 1

   references
   rnaseq
   chipseq

While the references workflow can be stand-alone, usually it is run as
a by-product of running the RNA-seq or ChIP-seq workflows. Here we will
focus on RNA-seq and ChIP-seq which share some common properties.

Where possible, we prefer to have rules use the normal command-line syntax for
tools (examples include rules calling samtools, deepTools bamCoverage, picard,
salmon).  However in some cases we use wrapper scripts. Situtations where we
use wrappers:

- Ensuring various aligners (HISAT2, Bowtie2, STAR, bwa) behave uniformly.
  These wrappers call the aligner, followed by samtools sort and view. The end
  result is that FASTQs go in, and a sorted BAM comes out.
- Tools with legacy dependencies like Python 2.7 that must be run in an
  independent environment (macs2, sicer, rseqc)
- R analyses (particularly spp and dupradar, which build up an R script
  incrementally before calling it).
- Tools that need complicated setup, or handling output files hard-coded by the
  tool (fastqc, fastq_screen).

In all cases, search for the string **NOTE:** in the Snakefile to read notes on
how to configure each rule, and make adjustments as necessary. You may see some
comments that say `# [TEST SETTINGS]`; you can ignore these, and see
:ref:`test-settings` for more info.

.. note:: 

    If you have two different RNA-seq experiments, from different species, they
    have to be run separately. However, if downstream analyses will use them both
    then you would like to keep them in the same project. In this case, you can copy
    the ``workflows/rnaseq`` directory to two other directories:

    .. code-block:: bash

        cp -r workflows/rnaseq workflows/genome1-rnaseq
        cp -r workflows/rnaseq workflows/genome2-rnaseq

    Now, downstream analyses can link to and utilize results from these individual
    folders, while the whole project can remain self-contained.

Integrative analysis workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The integrative analysis workflows take input from the primary workflows and
tie them together.

The integrative analysis workflows are described in :ref:`integrative`:

.. toctree::
   :maxdepth: 2

   integrative

Features common to workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this section, we will take a higher-level look at the features common to
the primary analysis workflows.

- There is some shared code across the multiple Snakefiles. For instance,
  The directory ``../..`` is added to Python's path. This way, the ``../../lib``
  module can be found, and we can use the various helper functions there. This is
  also simpler than providing a `setup.py` to install the helper functions.

- The config file is hard-coded to be `config/config.yaml`. This allows the config file to be
  in the `config` dir with other config files without having to be specified on
  the command line, while also affording the user flexibility. For instance, a custom
  config can be specified at the command-line, using  ``snakemake
  --configfile <path to other config file>``.

- The config file is loaded using ``common.load_config``. This function resolves
  various paths (especially the references config section) and checks to see
  if the config is well-formatted.

- To make it easier to work with the config, a `SeqConfig` object is created. It
  needs that parsed config file as well as the patterns file (see
  :ref:`patterns-and-targets` for more on this). The act of creating this object
  reads the sample table, fills in the patterns with sample names, creates
  a reference dictionary (see ``common.references_dict``) for easy access to
  reference files, and for ChIP-seq, also fills in the filenames for the
  configured peak-calling runs. This object, called ``c`` for convenience, can be
  accessed to get all sort of information -- ``c.sampletable``, ``c.config``,
  ``c.patterns``, ``c.targets``, and ``c.refdict`` are frequently used in rules
  throughout the Snakefiles.

- Various files can be used to specify cluster-specific parameters if the workflows
  are being run in a high-performance cluster environment. For example, a config file
  ``config/clusterconfig.yaml`` can be used to specify global and rule-specific
  memory and disk-space requirements for the Snakefile to use at run-time. For more
  details, see :ref:`cluster`.

Next Steps
~~~~~~~~~~

Next we look at :ref:`config` for details on how to configure specific workflows,
before going into the implemented workflows:

- Primary analysis workflows
   - :ref:`references`
   - :ref:`rnaseq`
   - :ref:`chipseq`

- :ref:`integrative`

