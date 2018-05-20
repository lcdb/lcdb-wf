.. _running-the-tests:

Running the tests
=================

This section describes setting up for and running the example data. This
reproduces the steps that are performed during the automated tests on Circle CI
(https://circleci.com/gh/lcdb/lcdb-wf).

The example run takes up about 360 MB of space and runs in about 15 mins on
2 cores.

The following will probably just be done once to verify the tests are working;
you won't need this when analyzing your own data.

.. note::

    The environment you created above needs to be activated for these steps.


1. Download example data
------------------------

This will download the example data from our `test data repository
<https://github.com/lcdb/lcdb-test-data>`_ into the directories
``workflows/{references,rnaseq,chipseq}/data``::

    python ci/get-data.py


.. _test-settings:

1a. A note about test settings
------------------------------

.. warning::

    Default configuration assumes a machine with large amounts of RAM. Running
    the workflows as-is on a single machine with limited RAM may cause all RAM
    to be consumed! Use ``run_test.sh`` as described below to avoid this.

A major benefit of ``lcdb-wf`` is that the code undergoes automated testing
(currently on `CircleCI <https://circleci.com/gh/lcdb>`_). However this test
environment only has 2 cores and 2GB RAM. To accommodate this but to also allow
the workflows to run in their entirety in a reasonable time frame, We developed
a small representative `test dataset <https://github.com/lcdb/lcdb-test-data>`_
from real-world data. We also need to make settings to the worflows, in
particular to set the Java VM memory to only 2GB for Java tools like Picard and
FastQC.

We had to make a design decision: should the "default" state of
the workflows reflect production-ready (high-RAM) settings, or reflect
test-ready (low RAM) settings? We chose to have the default to be real-world,
production-ready settings, because we want to minimize the editing (and
therefore possibility of introducing errors) in production.

So to run tests, we need to make some adjustments. In each workflows directory,
a ``run_test.sh`` script handles this. This script runs a preprocessor,
``ci/preprocessor.py``, which looks for specially-formatted comments in the
workflows, swaps out production settings for test settings, and writes the
results to stdout. In production, especially when running on a cluster, there's
no need to do this.

The ``run_test.sh`` simply passes all arguments on to Snakemake. Take a look at
the script to see what it's doing, and see the examples below for usage.

2. Try the RNA-seq workflow with example data
---------------------------------------------


With the `lcdb-wf` environment activated, change to the RNA-seq workflows
directory:

.. code-block:: bash

    cd workflows/rnaseq

First, run in dry-run mode which will print out the jobs to be run.  The
arguments will be described later, this is just to get things running:

.. code-block:: bash

    ./run_test.sh -n --use-conda --configfile ../../include/reference_configs/test.yaml

If all goes well, you will get lots of output ending with a summary of the
number of jobs that will be run. Then, use the same command but remove the
``-n``, and optionall include the ``-j`` argument to specify the number of
cores to use, for example ``-j 8`` if you have 8 cores on your machine (this
example just uses 2 cores):

.. code-block:: bash

    ./run_test.sh -j 2 --use-conda --configfile ../../include/reference_configs/test.yaml

This will take ~15 minutes to run.

Briefly, this workflow first imports the references workflow, which downloads
genome sequence and reference files and builds indexes as necessary (HISAT2
genome index, salmon transcriptome index, bowtie2 index for rRNA, GTF file of
gene annotations) and then carries on with the RNA-seq workflow.

The RNA-seq workflow includes the standard mapping, counting, and differential
expression stages, as well as many quality-control steps. See :ref:`rnaseq` for
more details.

After the workflow runs, here are some useful points of interest in the output:

    - ``data/rnaseq_samples/*``: sample-specific output. For example,
      individual BAMs and bigWig files can be found here
    - ``data/aggregation/multiqc.html``:  MultiQC report.
    - ``downstream/rnaseq.html``: Differential expression results generated
      from running the ``downstream/rnaseq.Rmd`` RMarkdown file.

See :ref:`rnaseq` for details.

Run the ChIP-seq workflow with example data
-------------------------------------------

With the `lcdb-wf` environment activated, from the top-level directory of the
repo, change to the ``workflows/chipseq`` directory:

.. code-block:: bash

    cd workflows/chipseq

First, run in dry-run mode which will print out the jobs to be run.  The
arguments will be described later, this is just to get things running:

.. code-block:: bash

    ./run_test.sh -n --use-conda --configfile ../../include/reference_configs/test.yaml

If all goes well, you will get lots of output ending with a summary of the
number of jobs that will be run. Then, use the same command but remove the
``-n``, and optionall include the ``-j`` argument to specify the number of
cores to use, for example ``-j 8`` if you have 8 cores on your machine (this
example just uses 2 cores):

.. code-block:: bash

    ./run_test.sh -j 2 --use-conda --configfile ../../include/reference_configs/test.yaml

Like the RNA-seq workflow, the ChIP-seq workflow includes the
``workflows/references/Snakemake`` workflow, so that genome fastas are
downloaded and indexes built as necessary, before continuing on to the ChIP-seq
workflow.

Points of interest:

    - ``data/chipseq_samples/*``: sample-specific output. Individual BAM files
      for a sample can be found here.
    - ``data/chipseq_merged/*``: technical replicates merged and re-deduped, or
      if only one tech rep, symlinked to the BAM in the samples directory
    - ``data/chipseq_peaks/*``: peak-caller output, including BED files of
      called peaks and bedGraph files of signal as output by each algorithm
    - ``data/chipseq_aggregation/multiqc.html``: MultiQC report

See :ref:`chipseq` for details.

Run the references workflow with example data
---------------------------------------------

This is optional; parts of this workflow were actually run automatically as
needed for the RNA-seq and ChIP-seq workflows. However, running this workflow
on its own can be useful for setting up a new site, as it will build all
configured references the config file provided to it (as opposed to only
building the references specifically requested by either the ChIP-seq or
RNA-seq workflows).

From the top-level of the repo, change to the ``workflows/references`` directory:
.. code-block:: bash

    cd workflows/references

First, run in dry-run mode which will print out the jobs to be run.  The
arguments will be described later, this is just to get things running:

.. code-block:: bash

    ./run_test.sh -n --use-conda --configfile ../../include/reference_configs/test.yaml

If all goes well, you will get lots of output ending with a summary of the
number of jobs that will be run. Then, use the same command but remove the
``-n``, and optionall include the ``-j`` argument to specify the number of
cores to use, for example ``-j 8`` if you have 8 cores on your machine (this
example just uses 2 cores):

.. code-block:: bash

    ./run_test.sh -j 2 --use-conda --configfile ../../include/reference_configs/test.yaml


See :ref:`references` for details.


Next steps
----------
See :ref:`config` for how to configure the workflows to work on your own data
and how to configure for your system.

See the :ref:`rnaseq`, :ref:`chipseq`, and :ref:`references` sections for more
details on the above workflows, and then the :ref:`external`, :ref:`figures`,
and :ref:`colocalization` sections for other workflows that can be used for
downstream analysis and integrating published data with newly-generated
results.
