.. _running-the-tests:

Testing the installation
========================
This section describes how to set up and run the example data.
It is useful for verifying everything is working correctly. This
reproduces the steps that are performed during the automated tests
on `Circle CI<https:/circleci.com>`_. You can see the latest test
results `here<https://circleci.com/gh/lcdb/lcdb-wf/tree/master>`_.

The example run takes up about 360 MB of space and runs in about 15 mins on
2 cores.

.. note::

   The ``deploy.py`` script specifically **excludes** the various test files,
   so the commands below must be run in a full clone of the repo, not in
   a directory in which lcdb-wf has been deployed.

Create conda envs
-----------------

This assumes you have set up the `bioconda channel
<https://bioconda.github.io>`_ properly.

.. code-block:: bash

   conda create -p ./env --file requirements-non-r.txt

.. code-block:: bash

   conda create -p ./env-r --file requirements-r.txt

We **highly recommend** using conda for isolating projects and for analysis
reproducibility. If you are unfamiliar with conda, we provide a more detailed look
at:

.. toctree::
   :maxdepth: 2

   conda


Activate the main env
---------------------

Depending on how you have set up conda, either

.. code-block:: bash

   conda activate ./env

or

.. code-block:: bash

   source activate ./env

Download example data
---------------------

This will download the example data from our `test data repository
<https://github.com/lcdb/lcdb-test-data>`_ into the directories
``workflows/{references,rnaseq,chipseq}/data``:

.. code-block:: bash

    python ci/get-data.py


.. _test-settings:

A note about test settings
--------------------------

.. warning::

    The default configuration assumes a machine with large amounts of RAM.
    Running the workflows as-is on a single machine with limited RAM may cause
    all RAM to be consumed! Use ``run_test.sh`` as described below to avoid
    this.

A major benefit of ``lcdb-wf`` is that the code undergoes automated testing on
`CircleCI <https://circleci.com/gh/lcdb>`_. However this test environment only
has 2 cores and 2GB RAM. To accommodate this, we developed a small
representative `test dataset <https://github.com/lcdb/lcdb-test-data>`_ from
real-world data.This allows the workflows to run in their entirety in a reasonable time frame.
We also need to adjust specific settings to the workflows, e.g.
we set the Java VM memory to only 2GB for Java tools like Picard and FastQC.

We had to make a design decision about the “default” state of the workflows:
should the workflows reflect production-ready (high-RAM) settings, or reflect
test-ready (low RAM) settings? We chose to have the default to be real-world,
production-ready settings, because we want to minimize the edits required
(and therefore possibility of introducing errors!) for running on real data.

What this all means is that if we want to run tests, we need to make some
adjustments. In each workflows directory, a ``run_test.sh`` script handles
this. This script runs a preprocessor, ``ci/preprocessor.py``, which looks for
specially-formatted comments in the workflows. It swaps out production settings
for test settings, and writes the results to a new ``Snakefile.test`` file that
is then run. In production, especially when running on a cluster, there's no
need to do this.


See the docstring in the ``ci/preprocessor.py`` for details on how this works.

The ``run_test.sh`` simply passes all arguments on to Snakemake. Take a look at
the script to see what it's doing, and see the examples below for usage.

Run the RNA-seq workflow with example data
------------------------------------------

With the `lcdb-wf` environment activated, change to the RNA-seq workflows
directory:

.. code-block:: bash

    cd workflows/rnaseq

First, run in dry-run mode which will print out the jobs to be run.  The
arguments will be described later, this is just to get things running:

.. code-block:: bash

    ./run_test.sh -n --use-conda

If all goes well, you will get lots of output ending with a summary of the
number of jobs that will be run. Then, use the same command but remove the
``-n``, and optionally include the ``-j`` argument to specify the number of
cores to use, for example ``-j 8`` if you have 8 cores on your machine (this
example just uses 2 cores):

.. code-block:: bash

    ./run_test.sh -j 2 --use-conda

This will take ~15 minutes to run.

Then activate the R environment (this assumes you're still in the
``workflows/rnaseq`` subdirectory):

.. code-block:: bash

    conda activate env-r   # or source activate env-r

You can either start an R interpreter and run:

.. code-block:: bash

    rmarkdown::render('downstream/rnaseq.Rmd')

or from the terminal:

.. code-block:: bash

    Rscript -e "rmarkdown::render('downstream/rnaseq.Rmd')"

After the workflow runs, here are some useful points of interest in the output:

    - ``data/rnaseq_samples/*``: sample-specific output. For example,
      individual BAMs and bigWig files can be found here
    - ``data/aggregation/multiqc.html``:  MultiQC report.
    - ``downstream/rnaseq.html``: Differential expression results generated
      from running the ``downstream/rnaseq.Rmd`` RMarkdown file.

See :ref:`rnaseq` and :ref:`config` for more details.

Run the ChIP-seq workflow with example data
-------------------------------------------

To run the ChIP-Seq workflow, follow the same steps as above but
with the workflow directory updated to ``workflows/chipseq``.
The most notable difference here is that the downstream analysis
in R (e.g. the ``rmarkdown::render`` step)  is not run.

Points of interest after running the ChIP-seq workflow:

    - ``data/chipseq_samples/*``: sample-specific output. Individual BAM files
      for a sample can be found here.
    - ``data/chipseq_merged/*``: technical replicates merged and re-deduped, or
      if only one tech rep, symlinked to the BAM in the samples directory
    - ``data/chipseq_peaks/*``: peak-caller output, including BED files of
      called peaks and bedGraph files of signal as output by each algorithm
    - ``data/chipseq_aggregation/multiqc.html``: MultiQC report

See :ref:`chipseq` for more details.


Exhaustive tests
----------------

The file ``.circleci/config.yml`` configures all of the tests that are run on
CircleCI. There's a lot of configuration happening there, but look for the
entries that have ``./run_test.sh`` in them to see the commands that are run.

Next steps
----------

Now that you have tested your installation of ``lcdb-wf`` you can learn about the
different workflows implemented here at the :ref:`workflows` page and see details
on configuration at :ref:`config`, before getting started on your analysis.

In addition, :ref:`setup-proj` explains the process of deploying ``lcdb-wf``
to a project directory.
