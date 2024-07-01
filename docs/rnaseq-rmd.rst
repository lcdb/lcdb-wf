.. _downstream-detailed:

Detailed documentation of RNA-Seq downstream
============================================

Here we describe in detail the downstream analysis of RNA-Seq data performed in :file:`workflows/rnaseq/downstream/rnaseq.Rmd`. 

This page has one section per named chunk. For example, if :file:`rnaseq.Rmd`
has the following code:

.. code-block:: r

    ```{r load_libraries}
    library(DESeq2)
    library(dplyr)
    ```

Then you can find the corresponding documentation on this page under the
``load_libraries`` heading.

The file :file:`ci/ensure_docs.py` double-checks to make sure all chunks are
documented and all documentation corresponds to a chunk, as part of the testing
framework.

.. _rnaseqrmd:


``global_options``
------------------
This chunk sets global rmarkdown options. Some of the lines provide
a mechanism for all cached chunks to have their cache invalidated if either of
the filenames' modification times have changed. For example, the following
argument to ``knitr::opts_chunk$set`` will inject the option
``cache.extra_file_dep_1`` into all chunks:

.. code-block:: r

    cache.extra_file_dep_1=file.info('../config/sampletable.tsv')$mtime,

Therefore if the sample table of the RNA-Seq workflow has been updated since
the last time the analysis script was run, the modification time (mtime) will
be changed, so the ``cache.extra_file_dep_1`` value will be different than it
was before for every chunk, and so every chunk will be re-run.

``lcdbwf``
----------

Loads the ``lcdbwf`` R package, stored in ``../../../lib/lcdbwf/R``. This chunk
is not cached and fully reloads the package each time using
``devtools::load_all``, so any changes to the code in that package will show up
when this file is rendered. Documentation is also automatically re-generated.

Note that throughout this RMarkdown file, functions from this package will use
the ``lcdbwf::`` namespace prefix to be explicit about where that function is
coming from.

``config``
----------

This chunk loads the config files :file:`config.yaml` and :file:`text.yaml`.
This chunk is not cached, so any changes in the config files automatically show
up here.

See the :file:`config.yaml` file for configuring the code.

See the :file:`text.yaml` file for editing the explanatory text.

This chunk also configures the parallelization options in a chunk that is not
cached.

``libraries``
-------------

Standard loading of the library dependencies used throughout the code.

``coldata_setup``
-----------------
This chunk loads the sample table.

Use this chunk to add additional columns to your sampletable.

The only requirement is that the rownames correspond to sample names.

``dds_initial``
---------------

Constructs the initial dds object and the variance stabilized counts.

Note that the design used is ``~1`` and the call to
``varianceStabilizingTransformation`` uses ``blind=TRUE``, so you don't need to
change anything here.

Also note that the entire ``config`` object is passed to the ``make_dds``
function, which reads options like whether to strip version numbers off of gene
IDs or whether (and how) to collapse technical replicates.

``print_coldata``
-----------------

Simply prints the colData for reference -- excluding columns that might be in
there that would clutter the output.

``sample_heatmap``
------------------

This chunk creates a clustered heatmap of sample distances. The columns
specified in the config file's "covariates_for_plots" item will show up as
colors along the right side.


``pca``
-------

Creates PCA plots, one tab per entry in the config file's
"covariates_for_plots" item. Each plot has the same points, but the colors
differ. This can help assess the experimental design and set expectations on
how many differentially expressed genes one may find.

These are interactive plots, and hovering over a point indicates the sample.


``sizefactors``
---------------

This chunk makes diagnostic plots. In general, we expect sizeFactors to
correlate with total read count. When it doesn't, it can indicate that a small
number of genes are very highly expressed.


.. _dds_list:

``dds_list``
------------

This chunk sets up the :term:`dds` objects to be used in the `results` section
below for differential expression detection.

You may need different ``dds`` objects for testing different models, or perhaps
removing outlier samples. If you have technical replicates you might need to
combine them, and you might need to remove gene version identifiers. You might
want to use salmon instead of featureCounts. These would need to be done for
each ``dds``, requiring code duplication.

After working on many complex and/or messy experimental designs, we have
settled on the approach of a named list of ``dds`` objects, where later code
refers to these objects by their name in the list.

**The** ``results`` **chunk below expects such a list.**

The simplest example is the following where we create a single ``dds`` and put
it into a list.

.. code-block:: r

   dds <- DESeqFromCombinedFeatureCounts(
      '../data/rnaseq_aggregation/featurecounts.txt',
      sampletable=colData,
      design=~group)
   dds <- DESeq(dds, parallel=parallel)

   dds.list <- list(main=dds)

Now imagine a case where we want to remove a replicate that we think is an
outlier, but we still want to compare it to the results when it is included.
Let's say we also need to collapse the technical replicates. Such code would
look like this:

.. code-block:: r

   # The long way...
   #
   # First object with all replicates
   dds1 <- DESeqFromCombinedFeatureCounts(
      '../data/rnaseq_aggregation/featurecounts.txt',
      sampletable=colData,
      design=~group)
   dds1 <- collapseReplicates(dds1, 'biorep')
   dds1 <- DESeq(dds1, parallel=parallel)

   # Similar to above, but remove replicate 4
   dds2 <- DESeqFromCombinedFeatureCounts(
      '../data/rnaseq_aggregation/featurecounts.txt',
      sampletable=colData %>% filter(replicate!='rep4'),
      design=~group,
      # need subset_counts=TRUE if we want to automatically
      # subset the featureCounts to match the filtered colData
      # we provided.
      subset_counts=TRUE
      )
   dds2 <- collapseReplicates(dds, 'biorep')
   dds2 <- DESeq(dds2, parallel=parallel)

Based on our experience, as we add more ``dds`` objects the code gets more
error-prone. So for more complex use-cases, we have a function
``lcdbwf::make_dds``.

Here is how the code above would look using this method:

.. code-block:: r

   lst <- list(

      main=list(sampletable=colData, design=~group),

      no.rep.4=list(
         sampletable=colData %>% filter(replicate!='rep4'),
         design=~group,
         subset.counts=TRUE))
   )

   dds_list <- map(lst, lcdbwf::make_dds, config=config, parallel=config$parallel$parallel)

That is, first we create a list of lists (``lst``), and then we used ``map()`` to apply
the ``make_dds`` function to all items in the list. The collapsing of
replicates and other dds-creation configuration like stripping dotted version
names is determined by the config object which is passed along.

See the help for ``lcdbwf::make_dds`` for more details.

This chunk becomes a dependency of all of the ``results`` chunks below.

``dds_diagnostics``
-------------------

If configured, this chunk will run the diagnostics on the dds objects and show
tabbed reports on each dds.

``results_*``
-------------

.. note::

  This is where most of the customization needs to happen for each project.

This is actually a series of chunks where the bulk of the differential
expression analysis takes place.

For simple cases, you probably just need one of these. But for complex
experimental designs where you end up doing lots of contrasts, it can get time
consuming to run them every time you change the RMarkdown file.

The end result of these chunks is a single list containing DESeq2 results
objects and associated metadata in (sub)lists. Each of these sublists has:


- ``res``, the results object
- ``dds``, the string name in ``names(dds.list)``
- ``label``, a "nice" label which is used for headings and other output
- additional optional arguments that are passed along to ``DESeq2::results()``
  and/or ``DESeq2::lfcshrink()``

To continue our example from above, we might want to run the same contrast on
all samples (the "main" dds) and after removing replicate 4 (the "no.rep.4"
dds). To illustrate how additional arguments are used, let's imagine we also
want to use `ashr` as the shrinkage method for the second contrast.

Use the ``lcdbwf::make_results()`` function for this. This function is
a loose wrapper around ``DESeq2::results()`` and ``DESeq2::lfcshrink()`` that
adds some extra convenience when working with lists of dds objects, including
the detection of parallelization as set up in the config object. See the help
for ``lcdbwf::make_results()`` for more details.

By default, if no test argument is specified in the parameters for
``lcdbwf::make_dds`` (examples 1-4, rnaseq.Rmd, lines 164-187), the Wald test is
performed. When ``lcdbwf::make_results`` processes a Wald test dds object, it
detects the Wald test and expects a ``contrast`` or ``coef`` argument to specify which
p-values and log2FoldChange values to report.

DESeq2 also supports the nBinomLRT (LRT). Example 5 (rnaseq.Rmd, line 189)
demonstrates how to create a dds object with LRT data. Since the LRT tests
the removal of one or more terms from the design formula, a single
log2FoldChange column doesn't reflect the test's complexity. DESeq2's results
object is optimized for the Wald test, and when storing LRT results, it
maintains consistency in datastructure by choosing a single pair-wise comparison for
log2FoldChange values. To avoid confusion, ***we set all log2FoldChange values to
0 for LRT results***.

For more details, see the DESeq2 documentation: `DESeq2 Likelihood Ratio Test <https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#i-ran-a-likelihood-ratio-test-but-results-only-gives-me-one-comparison>`_.

.. _rules:

To take advantage of this infrastructure, we put each of those contrasts into
its own chunk **according to the following rules**:

- the chunk name must start with ``results_``
- the chunk is cached
- the chunk depends on the ``'dds_list'`` chunk
- the variable name starts with ``contr_[index]_``, and the rest of the variable name
  will be used as the name in the list. index is a string, containing 1 or more alphanumeric
  characters, and will be used as a sorting index for contrasts when generating output files.
  Note that the index string cannot contain "_".

Our example would look like the following. Note that we're showing the chunks
here because that will be come meaningful in a moment. They are shown as
comments here just to get the syntax highlighting to look OK.

.. code-block:: r

    # ```{r results_01, dependson='dds_list', cache=TRUE}
    contr_1_ko.vs.wt <- lcdbwf::make_results(
      dds_name='main',
      label='Using all samples',
      contrast=c('genotype', 'KO', 'WT')
    )
    # ```

    # ```{r results_02, dependson='dds_list', cache=TRUE}
    contr_2a_no.rep.4 <- lcdbwf::make_results(
      dds_name='no.rep.4',
      label='Removing replicate 4 and using ashr for shrinkage',
      contrast=c('genotype', 'KO', 'WT'),
      type='ashr'
    )
    # ```


When combined with the ``assemble_variables`` chunk described below, this
allows us to:

- retain caching at the level of individual contrasts
- combine all results into a single list used in later code while still respecting dependencies
- reduce the bookkeeping overhead

In our experience this scales well with very complex experimental designs with
lots of contrasts. For more information on creating complex contrasts, see
:ref:`contrast`. For more information on how these results are collected, see :ref:`assemble_variables`.


.. _assemble_variables:

``assemble_variables``
----------------------

If we had put all ``results()`` calls into the same chunk and cached that, then
a change anywhere in that chunk would invalidate the cache which would cause
all results to be regenerated. With many contrasts, this can get quite
time-consuming. An alternative would be to put each ``results()`` call into its
own chunk. But then we would need to keep track of dependencies and ensure
those dependencies were specified in downstream chunks.

For example, if you add a new results chunk and cache it, but forget to add
that chunk as a dependency in a later chunk, that later chunk will be
inconsistent and may even be missing the new results. Keeping track of this can
be error-prone.

Our solution is to set up the contrasts according to the :ref:`rules described
above <rules>`. By following those rules, the following becomes possible:

- we can detect all chunks creating results by looking for ``results_`` in the
  chunk name and automatically inject these into dependencies of future chunks.
- we can detect all results objects created by looking for variables starting
  with ``contr_[index]_``
- we can assemble all results objects into a list, and name each item in the
  list according to its variable name (minus the ``contr_[index]_``).
- we can alter the order of contrasts by simply modifying the index string in a single chunk.
  For example, if we have three contrasts contr_1_ko.vs.wt, contr_2_no.rep.4, and contr_3_no.rep.3, 
  we can change the order of contrasts simply by modifying one index string (ex: change contr_3_no.rep.3 to
  contr_1a_no.rep.3).

The ``assemble_variables`` chunk does all of this. The end result of this chunk
is a list of lists that is used by functions in the `lcdbwf` R package for
downstream work. For more details, see :term:`res_list`.

For each contrast (that is, each entry in `res_list`) the below chunks will
automatically create a DE results section including:

- a tabbed section using the label as a header
- summary table
- MA plot
- counts plots of top 3 up- and down-regulated genes
- p-value distribution
- exported results tables with links

.. _contrast:

Specifying contrasts
^^^^^^^^^^^^^^^^^^^^

Contrasts can be specified in three different ways.

.. note::

   In these examples, "control" and "treatment" are factor levels in the
   "group" factor (which was in the :term:`colData`), and the :term:`dds`
   object was created with the design ``~group``:

1. A character vector to the `contrast` parameter.

   This should be a three element vector:

   - the name of a factor in the design formula
   - name of the numerator for the fold change
   - the name of the denominator for the fold change. E.g.,

   .. code-block:: r

      res <- results(dds, contrast=c('group', 'treatment', 'control')

   That is, **the control must be last**.

2. `name` parameter for ``results()`` function call or `coef` parameter for
   ``lfcShrink()`` call

   `name` or `coef` should be one of the values returned by
   ``resultsNames(dds)`` that corresponds to the precomputed results. E.g.

   .. code-block:: r

      resultsNames(dds)
      # [1] "Intercept"  "group_treatment_vs_control"

      res <- results(dds, name='group_treatment_vs_control')

3. A numeric contrast vector with one element for each element in the
   ``resultsNames()`` function call. This is useful for arbitrary comparisons
   in multi-factor designs with a grouping variable.

   .. code-block:: r

      resultsNames(dds)
      # [1] "Intercept"  "group_treatment_vs_control"

      res <- results(dds, contrast=c(0, 1))


The most general way to specify contrasts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most general way to specify contrasts is with a numeric vector (third
option above).

Here is a worked example, using a two-factor experiment.

`group` encodes all combinations of a two-factor experiment, so we construct
a sampletable that looks like the following (here, showing 2 replicates per
group):

.. code-block::

   sample   genotype   condition   group
   1        A          I           IA
   2        A          I           IA
   3        B          I           IB
   4        B          I           IB
   5        A          II          IIA
   6        A          II          IIA
   7        B          II          IIB
   8        B          II          IIB



We can make arbitrary comparisons by fitting an 'intercept-less' model, e.g.
``design=~group + 0``, and numeric contrast vectors:

.. code-block:: r

   dds <- DESeqDataSetFromCombinedFeatureCounts(
       '../data/rnaseq_aggregation/featurecounts.txt',
       sampletable=colData,
       # NOTE: the design is now different
       design=~group + 0
   )
   dds <- DESeq(dds)

Check ``resultsNames``:

.. code-block::

   resultsNames(dds)
   # [1] "groupIA"  "groupIB"  "groupIIA"  "groupIIB"

So any numeric vectors we provide must be 4 items long. Here is how we can make
various contrasts with this experimental design. In each example, the
coefficients are indicated above the resultsNames to make it easier to see.

To compare IA and IB (that is, the genotype effect only in condition I):

.. code-block:: r

   #     1          -1         0           0
   # "groupIA"  "groupIB"  "groupIIA"  "groupIIB"

   res <- results(dds, contrast=c(1, -1, 0, 0)


Effect of genotype B (that is, disregard information about condition):

.. code-block:: r

   #     1          -1         1           -1
   # "groupIA"  "groupIB"  "groupIIA"  "groupIIB"

   res <- results(dds, contrast=c(1, -1, 1, -1)


Effect of condition II (that is, disregard information about genotype):

.. code-block:: r

   #     1           1         -1          -1
   # "groupIA"  "groupIB"  "groupIIA"  "groupIIB"

   res <- results(dds, contrast=c(1, 1, -1, -1)



Interaction term, that is, (IA vs IB) vs (IIA vs IIB). This is effectively ``(IA
- IB) - (IIA - IIB)``, which in turn becomes ``IA - IB - IIA + IIB``:

.. code-block:: r

   #     1          -1         -1          1
   # "groupIA"  "groupIB"  "groupIIA"  "groupIIB"

   res <- results(dds, contrast=c(1, -1, -1, 1)


``summary``
-----------

This chunk prints a high-level overview of all the contrasts.


``reportresults``
-----------------

This is the section that creates multiple, tabbed outputs for each of the
contrasts in the :term:`res_list`.

If the config specifies results diagnostics
(``config$toggle$results_diagnostics`` is TRUE), then this chunk will also run
the diagnostics. You can select just the ones you want diagnostics on using the
``config$plotting$diagnostics_results_names`` config option.

``upsetplots``
--------------

This chunk produces UpSet plots comparing the contrasts.

``excel``
---------

This chunk outputs an Excel spreadsheet with one contrast per sheet. Normalized
counts for each sample, from the respective dds object used for the contrast,
are also included on each sheet.

``write_output``
----------------

TSVs for each contrast's results are written to disk.

``combined_rds``
----------------

A single object is written as an .Rds file. This can then be used for
downstream visualization or it can be used as input to the functional
enrichment RMarkdown document.

``sessioninfo``
---------------

The output of sessionInfo records the versions of packages used in the analysis.

Glossary
--------
.. glossary::

   colData
      The metadata describing the samples. This is originally defined in the
      sampletable for the entire lcdb-wf run, is imported into rnaseq.Rmd, and
      may be subsequently modified.

   dds
      DESeq data set object. Typically this is incrementally added to, as in
      the DESeq2 vignette.

   vsd
      The variance-stabilized transformed version of the counts. Used for PCA,
      clustered heatmaps, and gene patterns.

   res_list
      A list, with one item per contrast. Each of those items in turn is a list
      of objects that together compose the contrast (dds name, results object, and
      label). This list-of-lists, which we call `res_list` for short, is used
      by functions in the `lcdbwf` R package for more downstream work, like
      gene patterns, functional enrichment, and Shiny apps.
