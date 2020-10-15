Overview of workflows
=====================

.. note::

   These workflows **are intended to be edited and customized by the user**.

   See :ref:`setup-proj` for recommendations on setting up these workflows in
   your project directory.


Orientation of the directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following is a high-level overview of the files used by the workflows. See
:ref:`config` for more specific configuration information.

Starting at the top level of the repo:

::

    [1]  ├── ci/
    [2]  ├── docs/
    [3]  ├── include/
    [4]  ├── lib/
    [5]  ├── README.md
    [6]  ├── requirements-non-r.txt
    [7]  ├── requirements-r.txt
    [8]  ├── workflows/
    [9]  └── wrappers/

1. ``ci`` contains infrastructure for continuous integration testing. You don't
   have to worry about this stuff unless you're actively developing `lcdb-wf`.

2. ``docs/`` contains the source for documentation. You're reading it.

3. ``include/`` has miscellaneous files and scripts that can be used by all
   workflows. Of particular note is the ``WRAPPER_SLURM`` script (see
   :ref:`cluster` for more) and the ``reference_configs`` directory (see
   :ref:`references` and :ref:`config` for more).

4. ``lib/`` contains Python modules used by the workflows.

5. ``README.md`` contains top-level info.

6. ``requirements-non-r.txt`` contains the package dependencies needed to run the
   workflows, and is used to set up a conda environment.

7. ``requirements-r.txt`` contains the package dependencies for R and various
   Bioconductor packages used in downstream analysis. See :ref:`conda-envs` for the
   rationale for splitting these.

8. ``workflows/`` contains one directory for each workflow; see below for details.

9. ``wrappers/`` contains Snakemake `wrappers
   <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#wrappers>`_,
   which are scripts that can use their own independent environment. See
   :ref:`wrappers` for more.


Workflow directories
--------------------
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
**integrative analysis**. The primary analysis workflows are:

- :ref:`references`
- :ref:`rnaseq`
- :ref:`chipseq`

The primary analysis workflows are generally used for transforming raw data
(fastq files) into usable results. For RNA-seq, that's differentially-expressed
genes (along with comprehensive QC and analysis). For ChIP-seq, that's called
peaks or differentially bound chromatin regions.

The integrative analysis workflows are:

- :ref:`colocalization`
- :ref:`external`
- :ref:`figures`

The integrative analysis workflows take input from the primary workflows and
tie them together.

Each workflow is driven by a ``Snakefile`` and is configured by plain text
`YAML <https://en.wikipedia.org/wiki/YAML>`_ and `TSV
<https://en.wikipedia.org/wiki/Tab-separated_values>`_ format files (see
:ref:`config` for much more on this).  In this section, we will take
a higher-level look at the features common to the primary analysis workflows.

Features common to workflows
----------------------------
There is some shared code across the multiple Snakefiles:

The directory ``../..`` is added to Python's path. This way, the ``../../lib``
module can be found, and we can use the various helper functions there. This is
also simpler than providing a `setup.py` to install the helper functions.

The config file is hard-coded to be `config/config.yaml`, but if you need to
this can be overridden by Snakemake from the commandline, using ``snakemake
--configfile <path to other config file>``. This allows the config file to be
in the `config` dir with other config files without having to be specified on
the command line, while also affording the user flexibility.

The config file is loaded using ``common.load_config``. This function resolves
various paths (especially the references config section) and checks to see
if the config is well-formatted.

To make it easier to work with the config, a `SeqConfig` object is created. It
needs that parsed config file as well as the patterns file (see
:ref:`patterns-and-targets` for more on this). The act of creating this object
reads the sample table, fills in the patterns with sample names, creates
a reference dictionary (see ``common.references_dict``) for easy access to
reference files, and for ChIP-seq, also fills in the filenames for the
configured peak-calling runs. This object, called ``c`` for convenience, can be
accesto get all sort of information -- ``c.sampletable``, ``c.config``,
``c.patterns``, ``c.targets``, and ``c.refdict`` are frequently used in rules
throughout the Snakefiles.

Cluster-specific settings
-------------------------
See :ref:`cluster` for details.


Primary analysis workflows
--------------------------
While the references workflow can be run stand-alone, usually it is run as
a by-product of running the RNA-seq or ChIP-seq workflows. See
:ref:`references` for details; here we will focus on RNA-seq and ChIP-seq which
share common properties.

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

    You can copy entire directories and keep them separate. As an example,
    imagine you have two different RNA-seq experiments. They are from two
    different species, and so have to be run separately. But you would like to
    keep them in the same project because downstream analysis will use them
    both.  In this case, you can copy the ``workflows/rnaseq`` directory to two
    other directories:

    .. code-block:: bash

        cp -r workflows/rnaseq workflows/genome1-rnaseq
        cp -r workflows/rnaseq workflows/genome2-rnaseq
