Overview of workflows
=====================

The following is a high-level overview of the files used by the workflows. See
:ref:`config` for more specific configuration information.

Starting at the top level of the repo:

::

    [1]  ├── ci/
    [2]  ├── docs/
    [3]  ├── include/
    [4]  ├── lib/
    [5]  ├── README.md
    [6]  ├── requirements.txt
    [7]  ├── workflows/
    [8]  └── wrappers/

1. ``ci`` contains infrastructure for continuous integration testing. You don't
   have to worry about this stuff unless you're actively developing `lcdb-wf`.

2. ``docs/`` contains the source for documentation. You're reading it.

3. ``include/`` has miscellaneous files and scripts that can be used by all
   workflows. Of particular note is the ``WRAPPER_SLURM`` script (see
   :ref:`cluster` for more) and the ``reference_configs`` directory (see
   :ref:`references` and :ref:`config` for more).

4. ``lib/`` contains Python modules used by the workflows.

5. ``README.md`` contains top-level info.

6. ``requirements.txt`` contains the package dependencies needed to run the
   workflows, and is used to set up a conda environment.

7. ``workflows/`` contains one directory for each workflow; see below for details.

8. ``wrappers/`` contains Snakemake `wrappers
   <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#wrappers>`_,
   which are scripts that can use their own independent environment. See
   :ref:`wrappers` for more.


Workflow directories
--------------------
Each workflow lives in its own directory:

::

    [1] ├── references
    [2] ├── rnaseq
    [3] ├── chipseq
    [4] ├── colocalization
    [5] ├── external
    [6] └── figures

There
are two general classes of workflows, the primary analysis and the downstream
analysis. The primary analysis workflows are:

1. :ref:`references`
2. :ref:`rnaseq`
3. :ref:`chipseq`

and the downstream workflows are:

4. :ref:`colocalization`
5. :ref:`external`
6. :ref:`figures`

Each workflow is driven by a ``Snakefile`` and is configured by plain text
`YAML <https://en.wikipedia.org/wiki/YAML>`_ and `TSV
<https://en.wikipedia.org/wiki/Tab-separated_values>`_ format files (see
:ref:`config` for much more on this).  In this section, we will take
a higher-level look at the features common to the primary analysis workflows.

Primary analysis workflows
--------------------------
While the references workflow can be run stand-alone, but usually it is run as
a by-product of running the RNA-seq or ChIP-seq workflows. See
:ref:`references` for details; here we will focus on RNA-seq and ChIP-seq which
share common properties.

When setting up an analysis, the three most important files in each workflow
directory are:


1. ``config/sampletable.tsv``. This configures which samples to run, their
   locations on disk, and other optional metadata for your experiment.
2. ``config/config.yaml``. This configures reference genomes to use and the
   location of your sampletable. For convenience this is a hard-coded filename,
   but see :ref:`config` for how to modify.
3. ``Snakefile``. This contains all the rules to run.

When the Snakefile is run, it loads the config file, which in turn loads the
configured sampletable. This information is used to tell Snakemake what files
need to be created.

Where possible, we prefer to have rules use the normal command-line syntax for
tools (examples include rules calling samtools, deepTools bamCoverage, picard,
salmon).  However in some cases we use wrapper scripts. Situtations where we
use wrappers:

- Aligners (HISAT2, Bowtie2). These wrappers call the aligner, followed by
  samtools sort and view such that FASTQs go in, and sorted BAM comes out.
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


