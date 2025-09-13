
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

Setting up a cluster profile
---------------------------

To run on a cluster, you need a Snakemake profile for your specific cluster
environment.

Typically this is set up by the cluster admins, but you can find some existing
profiles for common submission systems in the `Snakemake-Profiles project
<https://github.com/snakemake-profiles/doc>`__.

The profile, in turn, may require additional snakemake plugins to be
installed. These would need to go into the main environment; you can use the
``--additional-main`` argument for :file:`deploy.py` if/when you know the
package in advance. Otherwise, install the required plugin into the main
environment.

Then set the ``LCDBWF_SNAKEMAKE_PROFILE`` and ``LCDBWF_SNAKEMAKE_PROFILE_V8``
env vars to point to the <8 and >8 profiles respectively.


For NIH's Biowulf cluster
-------------------------

On NIH's Biowulf HPC cluster, use the `Biowulf profile
<https://github.com/NIH-HPC/snakemake_profile>`_. For recent versions
of lcdb-wf, you need the ``snakemake8`` branch of the profile:

.. code-block:: bash

    git clone https://github.com/NIH-HPC/snakemake_profile ~/snakemake_profile_v8
    cd ~/snakemake_profile_v8
    git checkout snakemake8

    # Set the environment variable (add to your ~/.bashrc for persistence)
    export LCDBWF_SNAKEMAKE_PROFILE_V8=~/snakemake_profile_v8

Running on a cluster
-------------------

In general, we should avoid running long-running Snakemake workflows on the
login node of a cluster. lcdb-wf comes with a wrapper script,
``include/WRAPPER_SLURM``, that runs Snakemake as a batch job on Slurm.

To run a workflow on a Slurm cluster, from the workflow directory (e.g.,
``workflows/rnaseq``), and with the main environment activated, run:

.. code-block:: bash

    sbatch ../../include/WRAPPER_SLURM


The environment is inherited, and the version of Snakemake is detected so that
the appropriate profile can be used. Then Snakemake is run. Via the profile,
the main Snakemake process started by the initial WRAPPER_SLURM batch job will
automatically submit jobs to the cluster.

You can keep an eye on :file:`Snakefile.log` to watch progress (e.g. `tail -f
Snakefile.log`).

For other cluster environments, consult the Snakemake documentation on:
- `Cluster execution <https://snakemake.readthedocs.io/en/stable/executing/cluster.html>`_
- `Cloud execution <https://snakemake.readthedocs.io/en/stable/executing/cloud.html>`_
