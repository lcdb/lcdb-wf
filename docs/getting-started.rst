.. _getting-started:

Getting started
===============

The only starting requirement is an installation of conda with the `bioconda
<https://bioconda.github.io>`_ channel set up.

.. note::

    `lcdb-wf` is tested and heavily used on Linux.

    It is likely to work on macOS as long as all relevant conda packages are
    available for macOS -- though this is not tested.

    It will **not** work on Windows due to a general lack of support of Windows
    in bioinformatics tools.

We use conda to install a separate environment for each project. This allows
long-running projects to keep using old versions of software if needed while
allowing newer versions in more recent projects.

To see the full list of software installed into these environments,
see the ``requirements-*.txt`` files at https://github.com/lcdb/lcdb-wf/.

For a more detailed discussion of conda, see :ref:`conda-envs`.

.. _setup-proj:

Setting up a project
--------------------

The general steps to use lcdb-wf in a new project are:

1. Deploy: download an instance of lcdb-wf to your project directory
2. Configure: set up samples table for experiments and edit configuration file
3. Run: activate environment and run the Snakemake file either locally or on a cluster

.. _deploy:

1. Deploying lcdb-wf
--------------------

lcdb-wf is designed to be customized for each project, so setting up a new
project to use lcdb-wf means copying over the files you need. There is
a ``deploy.py`` script that does this for you.

``deploy.py`` sets up a project and removes additional files (like these docs
and testing files) that are not necessary for an actual project. The reason to use this script is so you end up with a cleaner project directory. 

This script also writes a file to the destination called
``.lcdb-wf-deployment.json``. It stores the timestamp and details about what 
commit was used to deploy it. This tracks provenance of the code, so you can always
figure out what lcdb-wf commit your deployment originally started from.

To deploy a copy, you first need to get a copy. Here are three ways to do that:

Clone a repo using git:

.. code-block:: bash

   git clone git@github.com:lcdb/lcdb-wf.git

   # optionally check out a particular branch
   # git checkout BRANCHNAME

Or update an existing repo:

.. code-block:: bash

   BRANCH=master  # optionally change branch
   git checkout $BRANCH
   git pull origin $BRANCH

No git? Download a zip file:

.. code-block:: bash

   BRANCH=master  # optionally change branch
   wget https://github.com/lcdb/lcdb-wf/archive/$BRANCH.zip
   unzip $BRANCH.zip -d lcdb-wf


.. note::

   If you want to run the tests then don't deploy just yet -- see
   :ref:`running-the-tests` for details, and then come back here to deploy for
   an actual project.

Now that you have a copy of the code, you can deploy it to your project
directory. Run the deploy script in the newly cloned/updated repo. For help,
use ``-h``:

.. code-block:: bash

   python deploy.py -h

For example, to deploy the RNA-seq workflow to the ``myproj`` directory and
build the conda environments:

.. code-block:: bash

   python deploy.py --flavor rnaseq --dest myproj --build-envs

Copying over the files is fast; building the conda environments may take a few
minutes.

See :ref:`conda-envs` for more details on these.

2. Configure
------------

This step takes the most effort and the first time you set up a project it
will take some time to understand the configuration system.

- see :ref:`multiple-experiments` for advice on how to handle multiple experiments that are intended to be analyzed together
- see :ref:`conda-envs` for details on conda environments
- see :ref:`sampletable` for how to write a sampletable, which includes where to find raw data and contains the associated metadata
- see :ref:`config-yaml` for configuring each workflow

3. Run
------

Activate the main environment, go to the workflow you want to run, and run the
following:

.. code-block:: bash

    snakemake --dryrun

If all goes well, this should print a list of jobs to be run.

You can run locally, but this is NOT recommended. To run locally, choose the
number of CPUs you want to use (here, 8).

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

You can typically run simultaneous workflows when they are in different directories; see
:ref:`workflows` for details.

Next steps
~~~~~~~~~~

Next, we give a brief overview of the file hierarchy of ``lcdb-wf`` in the
:ref:`guide` page.
