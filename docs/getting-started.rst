.. _getting-started:

Getting started
===============

**Prerequisites:** `conda <https://docs.conda.io/en/latest/>`_ with the `bioconda <https://bioconda.github.io>`_ channel. See :ref:`conda-envs` if you need help setting this up.

.. note::

    `lcdb-wf` is only supported on Linux.

.. _setup-proj:

Quick start
-----------

1. **Deploy:** Copy workflow files to your project directory
2. **Configure:** Set up sample tables and edit config files
3. **Run:** Activate environment and run Snakemake

.. _deploy:

1. Deploy
---------

`lcdb-wf` is deployed by copying relevant files (Snakefiles, configs,
infrastructure) to your project directory. This creates
a ``.lcdb-wf-deployment.json`` file tracking the commit used.

**Option A: Download deployment script** (quickest method)

.. code-block:: bash

    wget https://raw.githubusercontent.com/lcdb/lcdb-wf/master/deploy.py
    python deploy.py \
      --dest analysis/project \
      --staging /tmp/lcdb-wf-tmp \
      --branch master \
      --flavor rnaseq \
      --clone \
      --build-envs

**Option B: Clone repo first** (use this to run tests)

.. code-block:: bash

   git clone https://github.com/lcdb/lcdb-wf /tmp/lcdb-wf
   cd /tmp/lcdb-wf
   python deploy.py \
     --dest analysis/project \
     --flavor rnaseq \
     --build-envs

See :ref:`conda-envs` for conda environment details.

2. Configure
------------
For your workflow of interest (see :ref:`workflows`), change to that directory
(e.g., :file:`workflows/rnaseq`).

- Edit :file:`config/sampletable.tsv` to reflect your samples and additional
  metadata. The required columns depend on the respective workflow type. See :ref:`sampletable`.
- Edit :ref:`config/config.yaml`. This is also workflow-specific, but at least
  points to the reference files (genome fasta). See :ref:`config`.

3. Run
------

.. warning::

    Some jobs require substantial RAM (e.g., 20 GB for typical MarkDuplicates;
    64 GB for building a STAR index for a mammalian genome). For
    MarkDuplicates, the Java VM will try to allocate this much RAM before
    starting, and will immediately crash if not enough is available. The STAR
    index building will continue to consume RAM and the machine may become
    sluggish and eventually crash.

Activate the environment and navigate to your workflow:

.. code-block:: bash

    conda activate ./env
    cd workflows/rnaseq
    snakemake --dryrun

**Local execution** (not recommended for typical projects):

.. code-block:: bash

    snakemake --use-conda -j 8

**Cluster execution** (recommended):

Use a `Snakemake profile
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ for
your cluster. For NIH's Biowulf, you can use the `NIH-HPC snakemake profile
<https://github.com/NIH-HPC/snakemake_profile>`__, using the ``snakemake8``
branch, like this:

.. code-block:: bash

    # One-time setup:
    git clone https://github.com/NIH-HPC/snakemake_profile ~/snakemake_profile
    (cd ~/snakemake_profile && git checkout snakemake8)

    # add this to ~/.bashrc for persistence
    export LCDBWF_SNAKEMAKE_PROFILE=~/snakemake_profile

    # Submit job to the Slurm cluster
    sbatch ../../include/WRAPPER_SLURM

For other clusters, see Snakemake's `cluster execution
<https://snakemake.readthedocs.io/en/stable/executing/cluster.html>`_ and
`cloud execution
<https://snakemake.readthedocs.io/en/stable/executing/cloud.html>`_ docs. See
:ref:`cluster` and :ref:`workflows` for more details.

4. Downstream analyses
----------------------

For RNA-seq, there is a comprehensive set of downstream analyses to run, see
:ref:`rnaseq-downstream`.

5. Review output
----------------

In the workflow's directory (e.g., :file:`workflows/rnaseq`) there will be
a :file:`data` directory with the output. See the respective details for each
workflow at :ref:`workflows`.
