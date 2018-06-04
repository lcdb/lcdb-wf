.. _downstream-rnaseq:

Downstream RNA-seq in R
=======================
In a typical RNA-seq analysis, it is relatively straightforward to go from raw
reads to read counts in features to importing them into R. After that however,
expression analysis gets a bit more complicated and highly depends on the
design of the experiment. We had to make a design decision: do we expose
various configuration options for an RNA-seq analysis in order to keep the
entire workflow driven by a config file, or allow the full flexibility of R?

We chose to use the full power and flexibility of R. We provide a fairly
complete RMarkdown file that is **intended to be edited and customized per
experiment**. This file is ``workflows/rnaseq/downstream/rnaseq.Rmd``, and it
is run by the ``rnaseq`` rule in the RNA-seq workflow, which renders it through
Knitr to create an HTML file. The rule's inputs are the featureCounts output,
the sample table, and the Rmd itself, so if any of these change the Rmd will be
re-run.

Many of the code chunks have the ``cached=TRUE`` option to speed up
re-rendering and make iterative development quicker. When everything's in
a final state, you may want to delete the ``rnaseq_cache`` directory and
re-run.

Editing ``rnaseq.Rmd``
----------------------
Search for the string ``NOTE:`` in ``rnaseq.Rmd`` to identify locations in the
RMarkdown file that should be edited, or at least checked, prior to running.

As with many analyses in R, the work is highly iterative. You may want to
consider using RStudio when working on the analysis, and then run the Snakefile
when complete.  To ensure that RStudio is using the same packages as the
workflows, you should set the ``RSTUDIO_WHICH_R`` environment variable.

The easiest way to do this is to activate the conda environment you're using
for the analysis, then export the identified location of R to that variable:

.. code-block:: bash

    source activate lcdb-wf
    export RSTUDIO_WHICH_R=$(which R)

On MacOS, you may additionally need the following:

.. code-block:: bash

    launchctl setenv RSTUDIO_WHICH_R $RSTUDIO_WHICH_R

Then run RStudio, which should pick up the conda environment's version of R and
should already have packages like DESeq2 installed in the environment.

Points of interest
------------------
Standard QC
~~~~~~~~~~~
Clustered correlation heatmaps, PCA plots, and clustered heatmaps of the most
variable genes are created, much like the `Bioconductor RNA-seq workflow
<https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html>`_.

Use Salmon counts or featureCounts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Both Salmon quantification and featureCounts are imported into ``rnaseq.Rmd``,
and you can decide which to use for any particular analysis.

Easy support for multiple contrasts, and comparing them
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For complex experimental designs, populate the ``res.list`` list with lists
containing the results, the DESeq object, and a nice label. For example:

.. code-block:: r

    res.list[['mut1']] <- list(
        res=lfcShrink(dds, contrast=c('group', 'mutant1', 'control')),
        dds=dds,
        label='Effect of mutant 1'
    )
    res.list[['time']] <- list(
        res=lfcShrink(dds.with.different.design, contrast=c('time', 't1', 't0')),
        dds=dds.with.different.design,
        label='Effect of time in WT'
    )

For each entry in ``res.list``, a section will be built with a header
reflecting the label, a summary table, counts plots of the top 3 up- and top
3 down-regulated genes, an M-A plot, a p-value distribution histogram, and
links to the exported results. In addition if there were multiple contrasts
with detected DE genes, these will be compared in UpSet plots to identify the
genes that are shared as DE across contrasts.

Help docs
~~~~~~~~~
For new users, or when distributing the output to collaborators who might
not be familiar with the plots containined in the report, a background and
help section are included. This can be disabled by deleting the "helpdocs"
chunk at the end of the file.

Helper functions
~~~~~~~~~~~~~~~~
The file ``helpers.Rmd`` contains lots of helper functions that are
included in ``rnaseq.Rmd`` as a child document. This way the main document
doesn't get cluttered with all the functions while you're editing it, but
the code for those functions is still included in the HTML output.

GO analysis
~~~~~~~~~~~
For each of the entries in ``res.list`` (as described above), GO analysis will
be performed using clusterProfiler, using GO and KEGG enrichment. These will be
shown in bubble plots, which each contrast in a different color so that GO
enrichment can be compared between contrasts on the same plot.
