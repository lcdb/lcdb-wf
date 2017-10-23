Getting started
===============
The following steps will run the example workflows on a relatively small test
data set to ensure that everything runs correctly. On 2 cores it takes about 15
minutes; use more cores to speed it up.

Once everything runs correctly, you can modify the files for your own
experiment.

bioconda (one-time setup)
-------------------------
Follow the instructions for setting up `bioconda <https://bioconda.github.io>`_.
This includes installing `conda` and setting up the channels.

Add the `lcdb` channel (one-time setup)
---------------------------------------

This enables the `lcdb` conda channel so that additional dependencies not
included in bioconda can be installed::

    conda config --add channels lcdb

Create a new conda environment (one-time setup)
-----------------------------------------------

This creates a top-level environment with Snakemake and other requirements. It
should be activated any time you'll be working with `lcdb-wf`. Here we're using
the name "lcdb-wf", but you can use anything::

    conda create -n lcdb-wf -y python=3.5 --file requirements.txt

Then activate the environment::

    source activate lcdb-wf

Or deactivate::

    source deactivate

Download example data (one-time setup)
--------------------------------------

This will download the example data to the directory ``data/``::

    python ci/get-data.py


Run the references workflow
---------------------------

With the `lcdb-wf` environment activated::

    snakemake -prs references.snakefile --configfile config/test_config.yaml --use-conda -j1

Adjust the ``-j`` argument to match the number of CPUs to run jobs in parallel and speed up the workflow.

The references are configured in `config/test_config.yaml` and that file
indicates the output will be in `references_data`, so you can inspect that
folder.

.. image:: references.png


Run the rnaseq workflow
-----------------------

Similarly::

    snakemake -prs rnaseq.snakefile --configfile config/test_config.yaml --use-conda -j1

.. image:: rnaseq.png

Points of interest:

    - ``data/rnaseq_samples/*``: sample-specific output
    - ``data/aggregation/multiqc.html``:  MultiQC report
    - ``downstream/rnaseq.html``: Differential expression results

Run the chipseq workflow
------------------------

::

    snakemake -prs chipseq.snakefile --configfile config/test_chipseq.yaml --use-conda -j1

Points of interest:

    - ``data/chipseq_samples/*``: sample-specific output
    - ``data/chipseq_merged/*``: technical replicates merged, or if only one tech rep, symlinked
    - ``data/chipseq_peaks/*``: peak-caller output, including BED files of
      called peaks and bedGraph files of signal as output by each algorithm
    - ``data/chipseq_aggregation/multiqc.html``: MultiQC report

Next steps
----------
- Add your original FASTQ files to ``data/rnaseq_samples/`` directories, likely
  using symlinks
- Edit ``config/sampletable.tsv`` and ``config/config.yml`` to reflect your
  experimental design
- Edit ``rnaseq.snakefile`` to disable any steps you don't need
- After running the workflow, you can customize ``downstream/rnaseq.Rmd`` for
  your particular analysis.

The authoritative source on the config file is ``config/test_config.yaml``. This
YAML-formatted file contains comments on what each item is used for. For
a "production" config file, with references set up for multiple genomes, see
``config/config.yml``.
