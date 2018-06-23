.. _cluster:

Running on a cluster
--------------------
The example commands in :ref:`getting-started` describe running Snakemake
locally. For larger data sets, you'll want to run them on an HPC cluster.
Snakemake `supports arbitrary cluster commands
<http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html>`_,
making it easy to run these workflows on many different cluster environments.

The default configuration we provide is specific to the NIH Biowulf cluster.
To run a workflow on Biowulf, from the workflow directory (e.g.,
``workflows/rnaseq``, run the following command::

    sbatch ../../include/WRAPPER_SLURM

The ``WRAPPER_SLURM`` script submits the main Snakemake process on a separate
node to avoid any restrictions from running on the head node. That main
snakemake process then submits each rule separately to the cluster scheduler.
As configured in that script, we specify ``config/clusterconfig.yaml`` as
containing the rule-specific cluster arguments.

That script also contains this Snakemake arguments::

    --jobname "s.{rulename}.{jobid}.sh" \
    --cluster 'sbatch {cluster.prefix} --cpus-per-task={threads}  --output=logs/{rule}.o.%j --error=logs/{rule}.e.%j' \

This means that each job is named after the rule and job id (the ``--jobname``
arg) and the stdout and stderr go to files in ``logs`` and are named after the
rule, followed by a ``.o`` or ``.e``, followed by the cluster job ID (the
``--cluster`` arg).
