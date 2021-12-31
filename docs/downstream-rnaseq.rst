.. _downstream:

RNA-Seq downstream analysis
===========================

In a typical RNA-seq analysis, it is relatively straightforward to go from raw
reads to read counts in features to importing them into R. After that however,
expression analysis gets a bit more complicated and highly depends on the
design of the experiment.

We attempted to strike the balance between simplicity -- where as much
configuration as possible takes place via a config file -- and flexibility
where the R code can be modified as needed depending on the project.

This file is ``workflows/rnaseq/downstream/rnaseq.Rmd``. It uses a separate
conda environment that just has the R dependencies. It is rendered via
``knitr`` to create an HTML file. The inputs for the rule are the featureCounts
output, the sample table, the ``lib/lcdbwf`` R package, and the Rmd.

.. warning::

   This RMarkdown file is **intended to be edited and customized per experiment**.


How to use this code
~~~~~~~~~~~~~~~~~~~~

1. Activate the ``env-r`` conda environment (created as part of setting up the
   `lcdb-wf` deployment)

2. Edit the :file:`workflows/rnaseq/downstream/config.yaml` file. It is
   heavily commented and should be self-explanatory.

3. Customize the contrasts you want to run (see below for details on this)

4. From the :file:`workflows/rnaseq/downstream` directory, run
   ``rmarkdown::render("rnaseq.Rmd")`` to get :file:`rnaseq.html`

Here are some additional notes:

- Many of the code chunks have the ``cache=TRUE`` option to speed up
  re-rendering and make iterative development quicker. When everything's in
  a final state, you may want to delete the ``rnaseq_cache`` directory and
  re-run.

- Many of the cached code chunks also specify a config argument. These config
  items are taken from the :file:`config.yaml` file living alongside the
  :file:`rnaseq.Rmd`. If a cached chunk specifies a config option, and the
  value in the config file changes, the chunk will be re-run because its cache
  is invalidated.

- As with many analyses in R, the work is highly iterative. You may want to
  consider using an interactive interpreter, either via the command line or
  RStudio. To ensure that RStudio is using the same packages as the workflows,
  you should set the ``RSTUDIO_WHICH_R`` environment variable.

  The easiest way to do this is to activate the conda environment you're using
  for the analysis, then export the identified location of R to that variable:

  .. code-block:: bash

      source activate lcdb-wf
      export RSTUDIO_WHICH_R=$(which R)

  On MacOS, you may additionally need the following:

  .. code-block:: bash

      launchctl setenv RSTUDIO_WHICH_R $RSTUDIO_WHICH_R

  Then run RStudio, which should pick up the conda environment's version of R and
  which will already have packages like DESeq2 installed in the environment.

More details
~~~~~~~~~~~~

For more detailed documentation, see :ref:`downstream-detailed`.

.. toctree::
   :maxdepth: 2

   rnaseq-rmd
