.. _getting-started:

Getting started
===============

The main prerequisite for `lcdb-wf` is `conda
<https://docs.conda.io/en/latest/>_`, with the `bioconda
<https://bioconda.github.io>`_. channel set up and the `mamba
<https://github.com/mamba-org/mamba>`_ drop-in replacement for conda.

If this is new to you, please see :ref:`conda-envs`.

.. note::

    `lcdb-wf` is tested and heavily used on Linux.

    It is likely to work on macOS as long as all relevant conda packages are
    available for macOS -- though this is not tested.

    It will **not** work on Windows due to a general lack of support of Windows
    in bioinformatics tools.

.. _setup-proj:

Setting up a project
--------------------

The general steps to use lcdb-wf in a new project are:

1. **Deploy:** download and run ``deploy.py``
2. **Configure:** set up samples table for experiments and edit configuration file
3. **Run:** activate environment and run the Snakemake file either locally or on a cluster

.. _deploy:

1. Deploying lcdb-wf
--------------------

Unlike other tools you may have used, `lcdb-wf` is not actually installed per
se. Rather, it is "deployed" by copying over relevant files from the `lcdb-wf`
repository to your project directory. This includes Snakefiles, config files,
and other infrastructure required to run, and excludes files like these docs
and testing files that are not necessary for an actual project. The reason to
use this script is so you end up with a cleaner project directory. 

This script also writes a file to the destination called
``.lcdb-wf-deployment.json``. It stores the timestamp and details about what
commit was used to deploy it. This tracks provenance of the code, so you can
always figure out what lcdb-wf commit your deployment originally started from.

There are a few ways of doing this.

Option 1: Download and run the deployment script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that you will not be able to run tests with this method, but it is likely
the most conveient method.

.. code-block:: bash

    BRANCH=master  # optionally change branch
    wget https://raw.githubusercontent.com/lcdb/lcdb-wf/$BRANCH/deploy.py

Run ``python deploy.py -h`` to see help. Be sure to use the ``--staging`` and
``--branch=$BRANCH`` arguments when using this method, which will clone the
repository to a location of your choosing. Once you deploy you can remove it. For example:

.. code-block:: bash

    python deploy.py \
      --dest analysis/project \
      --staging /tmp/lcdb-wf-tmp \
      --branch master \
      --flavor rnaseq \
      --build-envs

    # You can clean up the cloned copy if you want:
    # rm -rf /tmp/lcdb-wf-tmp

Option 2: Clone repo manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone a repo using git and check out the branch. Use this method for running
tests):

.. code-block:: bash

   BRANCH=master  # optionally change branch
   git clone https://github.com:lcdb/lcdb-wf /tmp/lcdb-wf
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

You can run locally, but this is NOT recommended. To run locally, choose the
number of CPUs you want to use with the ``-j`` argument as is standard for
Snakemake.

.. warning::

    If you haven't made any changes to the Snakefiles, be aware that the
    default configuration needs a lot of RAM. For example, the MarkDuplicates
    runs set 20 GB RAM for Java, and that's for each job. Adjust the Snakefiles
    accordingly if you don't have enough RAM available (search for "Xmx" to
    find the Java args that set memory).

.. code-block:: bash

    # run locally (not recommended)
    snakemake --use-conda -j 8

The recommended way is to run on a cluster. On NIH's Biowulf cluster, the way
to do this is to submit the wrapper script as a batch job:

.. code-block:: bash

    sbatch ../../include/WRAPPER_SLURM

and then monitor the various jobs that will be submitted on your behalf. See
:ref:`cluster` for more details on this.

Other clusters will need different configuration, but everything is standard
Snakemake. The Snakemake documentation on `cluster execution
<https://snakemake.readthedocs.io/en/stable/executing/cluster.html>`_ and
`cloud execution
<https://snakemake.readthedocs.io/en/stable/executing/cloud.html>`_ can be
consulted for running on your particular system.

You can typically run simultaneous workflows when they are in different
directories; see :ref:`workflows` for details.

Next steps
~~~~~~~~~~

Next, we give a brief overview of the file hierarchy of ``lcdb-wf`` in the
:ref:`guide` page.
