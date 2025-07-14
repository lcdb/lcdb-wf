.. _getting-started:

Getting started
===============

The main prerequisite for `lcdb-wf` is `conda <https://docs.conda.io/en/latest/>_`, with the `bioconda <https://bioconda.github.io>`_. channel set up and the `mamba <https://github.com/mamba-org/mamba>`_ drop-in replacement for conda installed.

If this is new to you, please see :ref:`conda-envs`.

.. note::

    `lcdb-wf` is tested and heavily used on Linux. It is only supported on
    Linux.

.. _setup-proj:

Setting up a project
--------------------

The general steps to use lcdb-wf in a new project are:

1. **Deploy:** download and run ``deploy.py`` to copy files into a project directory
2. **Configure:** set up samples table for experiments and edit configuration file
3. **Run:** activate environment and run the Snakemake file either locally or on a cluster

See :ref:`cluster` for detailed information on running on HPC clusters.

.. _deploy:

1. Deploying lcdb-wf
--------------------
Using `lcdb-wf` starts with copying files to a project directory, or
"deploying".

Unlike other tools you may have used, `lcdb-wf` is not actually installed per
se. Rather, it is "deployed" by copying over relevant files from the `lcdb-wf`
repository to your project directory. This includes Snakefiles, config files,
and other infrastructure required to run, and excludes files like these docs
and testing files that are not necessary for an actual project. The reason is
to use this script is so you end up with a cleaner project directory, compared
to cloning the repo directly.

This script also writes a file to the destination called
``.lcdb-wf-deployment.json``. It stores the timestamp and details about what
commit was used to deploy it. This tracks provenance of the code, so you can
always figure out what lcdb-wf commit your deployment originally started from.

There are a few ways of doing this.

Option 1: Download and run the deployment script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the most convenient method, although it does not allow running tests
locally.

.. code-block:: bash

    BRANCH=master  # optionally change branch
    wget https://raw.githubusercontent.com/lcdb/lcdb-wf/$BRANCH/deploy.py

Run ``python deploy.py -h`` to see help. Be sure to use the ``--staging`` and
``--branch=$BRANCH`` arguments when using this method, which will clone the
repository to a location of your choosing. Once you deploy you can remove the
script. For example:

.. code-block:: bash

    # Install any additional plugins required for running on your cluster
    EXTRA="snakemake-executor-plugin-cluster-generic"

    python deploy.py \
      --dest analysis/project \
      --staging /tmp/lcdb-wf-tmp \
      --branch $BRANCH \
      --flavor rnaseq \
      --clone \
      --build-envs \
      --additional_main=$EXTRA

    # You can clean up the cloned copy if you want:
    # rm -rf /tmp/lcdb-wf-tmp
    # and the downloaded script:
    # rm deploy.py


This will clone the full git repo to ``/tmp/lcdb-wf-tmp``, check out the master
branch (or whatever branch ``$BRANCH`` is set to), copy the files required for
an RNA-seq project over to ``analysis/project``, build the main conda
environment and the R environment, save the ``.lcdb-wf-deployment.json`` file
there, and then delete the temporary repo.

Option 2: Clone repo manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone a repo using git and check out the branch. Use this method for running
tests):

.. code-block:: bash

   BRANCH=master  # optionally change branch
   git clone https://github.com/lcdb/lcdb-wf /tmp/lcdb-wf
   cd /tmp/lcdb-wf
   git checkout $BRANCH

Then run the deploy script, ``python deploy.py -h`` to see usage info. Here is
an example for RNA-seq:

.. code-block:: bash

    python deploy.py \
      --dest analysis/project \
      --flavor rnaseq \
      --build-envs

.. note::

   If you want to run the tests then don't deploy just yet -- see
   :ref:`running-the-tests` for details, and then come back here to deploy for
   an actual project.


.. note::

    See :ref:`conda-envs` for more details on the conda environment building.

2. Configure
------------

This step takes the most effort. The first time you set up a project it
will take some time to understand the configuration system.

- see :ref:`sampletable` for how to write a sampletable, which includes where to find raw data and contains the associated metadata
- see :ref:`config-yaml` for configuring each workflow
- see :ref:`multiple-experiments` for advice on how to handle multiple experiments that are intended to be analyzed together

3. Run
------

Activate the main environment and go to the workflow you want to run. For
example if you have deployed and configured an RNA-seq run, then do:

.. code-block:: bash

    conda activate ./env
    cd workflows/rnaseq

and run the following:

.. code-block:: bash

    snakemake --dryrun

If all goes well, this should print a list of jobs to be run.

You can run locally, but this is NOT recommended for a typicaly RNA-seq
project. To run locally, choose the number of CPUs you want to use with the
``-j`` argument as is standard for Snakemake.

.. warning::

    If you haven't made any changes to the Snakefiles, be aware that the
    default configuration needs a lot of RAM. For example, the MarkDuplicates
    runs set 20 GB RAM for Java, and that's for each job. Adjust the Snakefiles
    accordingly if you don't have enough RAM available (search for "Xmx" to
    find the Java args that set memory).

.. code-block:: bash

    # run locally (not recommended)
    snakemake --use-conda -j 8

The recommended way is to run on a cluster.

Running on a cluster requires a `Snakemake profile
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ that translates 
resource requirements into arguments for your cluster's batch system.

For example, on NIH's Biowulf cluster:

.. code-block:: bash

    sbatch ../../include/WRAPPER_SLURM

This submits Snakemake as a batch job, which then submits individual workflow jobs to the cluster.

See :ref:`cluster` for detailed setup instructions for different cluster environments, including:

- Setting up Snakemake profiles
- Installing required plugins
- Configuring environment variables
- Running on specific clusters like NIH's Biowulf

The :ref:`cluster` section also links to Snakemake's documentation for various execution environments.

You can typically run simultaneous workflows when they are in different
directories; see :ref:`workflows` for details.

Next steps
~~~~~~~~~~

Next, we give a brief overview of the file hierarchy of ``lcdb-wf`` in the
:ref:`guide` page.
