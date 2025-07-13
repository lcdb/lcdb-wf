
.. _config:


Configuration
=============

General configuration
~~~~~~~~~~~~~~~~~~~~~

The majority of the work in setting up a new project is in the configuration --
which samples to run, where the data files are located, which references are
needed, etc.

**The entry point for configuration** is in the ``config/config.yaml`` file found
in each workflow directory. See :ref:`config-yaml` for more.

.. toctree::
   :maxdepth: 2

   config-yaml

The **references section** of the config file configures the genomes,
transcriptomes, and annotations to be used. See :ref:`references-config` for more.

.. toctree::
   :maxdepth: 2

   references-config

The **sample table**, lists sample IDs, filenames, and other metadata. Its path
is specified in the config file. See :ref:`sampletable` for more.

.. toctree::
   :maxdepth: 2

   sampletable

A **patterns file** only needs to be edited
if you're doing custom work. It determines the patterns of files that will be
created by the workflow. See :ref:`patterns-and-targets` for more.

.. toctree::
   :maxdepth: 2

   patterns-targets

.. _cluster:

Running on a cluster
~~~~~~~~~~~~~~~~~~~~
The example commands in :ref:`getting-started` describe running Snakemake
locally. For larger data sets, you'll want to run them on an HPC cluster.
Snakemake `supports arbitrary cluster commands
<http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html>`_,
making it easy to run these workflows on many different cluster environments.

Snakemake and these workflows are designed to decouple the code from the
configuration. Each rule has resources specified. When running with
a cluster-specific `Snakemake profile
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_,
these resources are translated into cluster-specific commands.

For example, if runnng NIH's Biowulf HPC cluster, use the `Biowulf profile
<https://github.com/NIH-HPC/snakemake_profile>`_. When using a recent version
of lcdb-wf, you will need the ``snakemake8`` branch of this repo (since pre-v8
snakemake profiles are incompatible).

Generally, you shouldn't run long-running tasks on a login node of a cluster,
and this includes long-running Snakemake workflows. So lcdb-wf comes with
a wrapper script, ``include/WRAPPER_SLURM``, that runs Snakemake which can be
submitted to a compute node on a Slurm cluster.

For example, to run a workflow on a Slurm cluster, from the workflow directory
(e.g., ``workflows/rnaseq``, run the following command::

    sbatch ../../include/WRAPPER_SLURM

The ``WRAPPER_SLURM`` script submits the main Snakemake process on a separate
node to avoid any restrictions from running on the head node. That main
Snakemake process then submits each rule separately to the cluster scheduler.
