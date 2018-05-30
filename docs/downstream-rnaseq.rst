.. _downstream-rnaseq:

Downstream RNA-seq in R
=======================
In a typical RNA-seq analysis, going from raw reads to read counts in features
and importing them into R is relatively straightforward. After that, expression
analysis gets a bit more complicated, depending on the design of the
experiment.

Rather than expose configuration options in a config file, we instead support
the full power and flexibility of R by providing a fairly complete RMarkdown
file that is intended to be edited and customized per experiment. This file is
``workflows/rnaseq/downstream/rnaseq.Rmd``, and it is run by the ``rnaseq``
rule in the RNA-seq workflow by rendering it through Knitr to create an HTML
file. The rule's inputs are the featureCounts output, the sample table, and the
Rmd itself, so if any of these change the Rmd will be re-run.

Many of the code chunks have the ``cached=TRUE`` option to speed up
re-rendering and make iterative development quicker. When everything's in
a final state, you may want to delete the ``rnaseq_cache`` directory and
re-run.

Configuration
-------------
Similar to the Snakefiles, search for the string ``NOTE:`` in ``rnaseq.Rmd``
to identify locations in the RMarkdown file that should be edited, or at least
checked, prior to running.

As with many analyses in R, the work is highly iterative. You may want to
consider using RStudio when working on the analysis, and then run the Snakefile
when complete.  To ensure that RStudio is using the same packages as the
workflows, you should set the ``RSTUDIO_WHICH_R`` environment variable.

The easiest way to do this is to activate the conda environment you're using
for the analysis, then export the identified location of R to that variable::

    source activate lcdb-wf
    export RSTUDIO_WHICH_R=$(which R)

On MacOS, you may additionally need the following::

    launchctl setenv RSTUDIO_WHICH_R $RSTUDIO_WHICH_R

Then run RStudio, which should pick up the conda environment's version of R and
should already have packages like DESeq2 installed.

Points of interest
------------------
For new users, or when distributing the output to collaborators who might not
be familiar with the plots containined in the report, a background and help
section are included.

The file ``helpers.Rmd`` contains lots of helper functions that are included in
``rnaseq.Rmd`` as a child document. This way the main document doesn't get
cluttered with all the functions while you're editing it, but the code for
those functions is still included in the HTML output.

