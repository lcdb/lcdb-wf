.. _cluster:

Running on a cluster
--------------------
The example commands in :ref:`getting-started` describe running Snakemake
locally. For larger data sets, you'll want to run them on an HPC cluster.
Snakemake `supports arbitrary cluster commands
<http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html>`_,
making it easy to run these workflows on many different cluster environments.

Snakemake, and these workflows, are designed to decouple the code from the
configuration. If you are running the workflows on NIH's Biowulf cluster, you
don't need to change anything. If you are running on a different cluster, you
should inspect the following files:

- `include/WRAPPER_SLURM`
- the `config/clusterconfig.yaml` files in each workflow directory you will be
  using (see :ref:`clusterconfig`)
- `lib/cluster_specific.py`. This module currently has a single function that,
  when called, will inspect the current environment variables and make any
  necessary changes, returning the temp dir. Other cluster-specific code may go
  here (see `cluster_specific`)

The default configuration we provide is specific to the NIH Biowulf cluster.
To run a workflow on Biowulf, from the workflow directory (e.g.,
``workflows/rnaseq``, run the following command::

    sbatch ../../include/WRAPPER_SLURM

The ``WRAPPER_SLURM`` script submits the main Snakemake process on a separate
node to avoid any restrictions from running on the head node. That main
snakemake process then submits each rule separately to the cluster scheduler.
As configured in that script, we specify ``config/clusterconfig.yaml`` as
containing the rule-specific cluster arguments.

That script also contains the Snakemake arguments::

    --jobname "s.{rulename}.{jobid}.sh" \
    --cluster 'sbatch {cluster.prefix} --cpus-per-task={threads}  --output=logs/{rule}.o.%j --error=logs/{rule}.e.%j' \

This means that each job will be named after the rule and job id (the
``--jobname`` arg) and the stdout and stderr go to files in ``logs`` and are
named after the rule, followed by a ``.o`` or ``.e``, followed by the cluster
job ID (the ``--cluster`` arg).

.. _cluster_specific:

TMPDIR handling
~~~~~~~~~~~~~~~
The top of each snakefile sets up a shell prefix that exports the TMPDIR
variable. The reason for this is that the NIH Biowulf cluster supports nodes
with temporary local storage in a directory named after the SLURM job ID. This
ID is not known ahead of time, but is stored in the ``SLURM_JOBID`` env var.

Since each rule executed on a cluster node calls the snakefile (see the job
scripts created by snakemake for more on this), we can look for the job ID and
set the tempdir appropriately. Upon setting ``$TMPDIR``, the Python
``tempfile`` module will use that directory to store temp files. Any wrappers
can additionally use ``$TMPDIR`` in shell commands and it will use this
directory.

Note that the default behavior -- if the ``SLURM_JOBID`` env var is not set --
is to set ``$TMPDIR`` to the default temp directory as documented in Python's
`tempfile module
<https://docs.python.org/3/library/tempfile.html#tempfile.gettempdir>`_.
However if you use these workflows on a different cluster, you may need to
provide a different function to return the job-specific temp directory.

.. _clusterconfig:

``clusterconfig.yaml`` files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `snakemake cluster configuration docs
<https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration>`_
describe how to configure cluster-specific settings so that rules will run on
the correct size node when submitting batch jobs. This is not needed if you are
running the workflows locally, but you may need to edit the
`config/clusterconfig.yaml` file found in each workflow directory.
