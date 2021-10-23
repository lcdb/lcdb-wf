.. _downstream-detailed:

Detailed documentation of RNA-Seq downstream
============================================

Here we describe in detail the downstream analysis of RNA-Seq data done using RMarkdown.
The code is broken into chunks or modules and the documentation follows the same
structure. For instance, if :file:`rnaseq.Rmd` has the following code:

.. code-block:: r

    ```{r load_libraries}
    library(DESeq2)
    library(dplyr)
    ```

Then you can find the corresponding documentation on this page under the
``load_libraries`` heading. Note that an RMarkdown code chunk must be named
in order to be documented. RMarkdown requires uniquely-named chunks, so we
can use them as uniquely-named headings.

The file :file:`ci/ensure_docs.py` double-checks to make sure all chunks are
documented and all documentation corresponds to a chunk, as part of the testing
framework.

To further modularize the code and because certain parts of the analysis are self-contained,
the downstream analysis code is broken up into three main components:

- `Primary analysis script`_: ``rnaseq.Rmd``
- `Gene patterns analysis`_: ``gene-patterns.Rmd``
- `Functional enrichment analysis`_: ``functional-enrichment.Rmd``

The ``gene_patterns.Rmd`` and ``functional-enrichment.Rmd`` scripts are called
from within the primary analysis Rmd in the genepatterns_ and functionalenrichment_
chunks below.

.. _rnaseqrmd:

Primary analysis script
~~~~~~~~~~~~~~~~~~~~~~~

``global_options``
------------------
This chunk sets global rmarkdown options. The following lines provide
a mechanism for all cached chunks to have their cache invalidated if either of
the filenames' modification times have changed.

.. code-block:: r

    cache.extra_file_dep_1=file.info('../config/sampletable.tsv')$mtime,
    cache.extra_file_dep_2 = file.info('../data/rnaseq_aggregation/featurecounts.txt')$mtime

In other words, if the sample table of the RNA-Seq workflow or the counts table from
the samples have been updated since the last time the analysis script was run, everything
will be rerun from scratch.

``load_helpers``
----------------

Helper functions are a set of useful utility functions used throughout the
``rnaseq.Rmd`` script that are stored in a package maintained in this
repository, in the ``lib/lcdbwf`` R package. Here we include these as a child
Rmd, so that the main document doesn't get cluttered with the functions, but
the code is still included in the HTML output.

This chunk refreshes documentation and loads (or re-loads, if you run it later)
the package using devtools.

``annotationhub_setup``
-----------------------

We use AnnotationHub for downloading annotations on the fly rather than
specifying OrgDbs in the requirements.txt of the conda environment. This allows
as much of our code as possible remain organism-agnostic. The package default
cache location is in the user's home directory. However this prevents other
people using the same environment from using the cached download. Therefore
here we reset the cache directory to
:file:`../../../include/AnnotationHubCache`.

+------------------------------+----------------------------------------------------------------------------------------------------------------------+
| var                          | description                                                                                                          |
+==============================+======================================================================================================================+
| ``annotation_genus_species`` | change this to what is relevant for this experiment. will search for and use the latest annotations for that species |
+------------------------------+----------------------------------------------------------------------------------------------------------------------+
| ``annotation_key_override``  | will use a specific key, if you know ahead of time what that is                                                      |
+------------------------------+----------------------------------------------------------------------------------------------------------------------+

``coldata_setup``
-----------------

This chunk prepares the :term:`colData` object that contains the metadata to be used
for creating the dds objects. The majority of the metadata is read in from the
sampletable. This chunk also handles things like converting to factors and
setting up the paths to Salmon output files.

+---------------------------+------------------------------------------------------------------------------------------------+
| var                       | description                                                                                    |
+===========================+================================================================================================+
| ``sample.table.filename`` | path to sampletable, generally you don't need to change this                                   |
+---------------------------+------------------------------------------------------------------------------------------------+
| ``strip.dotted.version``  | if TRUE, then Ensembl versions (the ".1" in "ENSG000012345.1") will be removed from gene names |
+---------------------------+------------------------------------------------------------------------------------------------+
| ``exclude.for.printing``  | default columns in the default sampletable that shouldn't necessarily be printed in tables     |
+---------------------------+------------------------------------------------------------------------------------------------+
| ``factor.columns``        | columns to ensure are factor types                                                             |
+---------------------------+------------------------------------------------------------------------------------------------+
| ``salmon.path.func``      | given a samplename ``x``, this function returns the sample's ``quant.sf`` file.                |
+---------------------------+------------------------------------------------------------------------------------------------+

This is a good place to put any modifications to the sample table (like factors
derived other columns).

If you are using this Rmd outside the context of lcdb-wf, you will need to
change the salmon output path patterns and the sampletable location.

.. topic:: Note on factors
   
   For the test data, "control" is the base level for the "group" factor. You will
   need to edit this as appropriate for your experimental design.


``salmon``
----------

If you don't want to use Salmon TPM, disable this chunk with ``eval=FALSE`` or
delete it entirely (and do the same with the next chunk).

``ddstxi``
----------

``design`` will likely need to be changed depending on your experimental
design.

This chunk creates separate ``dds.txi`` and ``vsd.txi`` objects to
differentiate them from the ones with no ``.txi`` that are created using
featureCounts.

Note we're using VST rather than rlog because the DESeq2 docs say they are
largely equivalent, and vst is substantially faster. Also note that since this
is exploratory analysis, we use ``blind=TRUE`` to ignore the design.

``dds_initial``
---------------
This initial :term:`dds` object will be used for exploratory data analysis, NOT
for differential expression. So the ``design`` should be something generic like
"group" even for complex experimental designs.

This chunk creates the initial :term:`dds` and :term:`vsd` objects that will be
used for exploratory data analysis.

``sample_heatmap``
------------------

This chunk creates a clustered heatmap of sample distances.

It can be helpful to add colors along the side to indicate different aspects of the
sample metadata. Any number of columns from the :term:`colData` can be provided
as ``cols.for.grouping``.

``pca``
-------

Create PCA plots, colored by possibly many different :term:`colData` columns
(specified using the ``groups`` list).

Each of the values in ``groups`` will have a corresponding interactive PCA plot
in a separate tab. This makes it easy to click through tabs to get a feel for
the structure of the data, and allows for hoving over a point to see the
metadata.

Note that plotting interactive plotly figures in a loop is not quite possible
(due to technical limitations) and so we have to use a workaround. Currently,
this workaround is to "manually" step through the loop, setting ``i`` to
a different integer and copy/pasting the same code multiple times.

``sizefactors``
---------------

To more easily investigate any outliers in these plots, you can optionally
attach columns from ``colData`` before plotting the scatterplot, e.g.:

.. code-block:: r

   color_by <- 'group'
   group_names <- tibble(name=dds$samplename, group=dds[[color_by]])
   trc_vs_sf <- full_join(sf, trc, by='name')

``parallel_config``
-------------------

By default we do not run in parallel, however this can be very useful in
experiments with many samples and complex designs. To run in parallel, manually
configure the parallel workers, set the number of cores, and set parallel to
TRUE:

.. code-block:: r

   parallel <- TRUE
   register(MulticoreParam(4))

Calls to ``DESeq()`` below will provide the argument ``parallel=parallel`` so no
other changes should be needed.

.. _dds_list:

``dds_list``
------------

This chunk sets up the :term:`dds` objects to be used in the `results` section
below for differential expression detection.

For simple cases, you can just create a :term:`dds` object and store it in
a single-item list. However this document is designed to work with quite
complex experimental designs, and we provide tools for conveniently working with
such complexity while hopefully reducing the possibility of errors.

You may need different ``dds`` objects for testing different models, or perhaps
removing outlier samples. If you have technical replicates you might need to
combine them, and you might need to remove gene version identifiers. You might
want to use salmon instead of featureCounts. These would need to be done for
each ``dds``, requiring code duplication.

After working on many complex and/or messy experimental designs, we have
settled on the approach of a named list of ``dds`` objects.

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

Now imagine we want to try removing a replicate that we think is an outlier, but
we still want to compare it to the results when using the full set of
replicates. Let's say we also need to collapse the technical replicates. Such
code would look like this:

.. code-block:: r

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
      # need subset.counts=TRUE if we want to automatically
      # subset the featureCounts to match the filtered colData
      # we provided.
      subset.counts=TRUE
      )
   dds2 <- collapseReplicates(dds, 'biorep')
   dds2 <- DESeq(dds2, parallel=parallel)

Based on our experience, as we add more ``dds`` objects the code gets more
error-prone. So for more complex use-cases, we have a function
``lcdbwf::make.dds``.

To use it, first we create a list of lists. The names of this list are useful
names you give each :term:`dds` object. Here, it's ``main`` and ``no.rep.4``.
For each of those names, the correspdonding values are lists with at least the
names ``sampletable`` and ``design`` which will be used to generate each
:term:`dds`. Aditional arguments to pass to ``DESeqFromCombinedFeatureCounts``,
like ``subset.counts=TRUE``, are provided in a separate ``args`` entry in the
list.

Then, we apply the ``make.dds`` function over that list:

.. code-block:: r

    map(lst, make.dds)

When doing so, we can optionally apply other arguments to every :term:`dds`
object in there. In the example below, we combine technical replicates on
biorep for every :term:`dds`, and use the same parallel argument for all of
them.

So the above example becomes the following:

.. code-block:: r

   lst <- list(
      main=list(sampletable=colData, design=~group),
      no.rep.4=list(
         sampletable=colData %>% filter(replicate!='rep4'),
         design=~group,
         args=list(subset.counts=TRUE))
   )

   dds.list <- map(lst, make.dds, combine.by='biorep', parallel=parallel)

Note the following:

- the file is set by default to be
  :file:`../data/rnaseq_aggregation/featurecounts.txt`. Use a different file on
  a dds-specific basis by including ``file="path/to/file.txt"``.
- we can supply additional args, like ``subset.counts=TRUE``, on a per-dds
  basis. If the sampletable is filtered, by default ``make.dds`` takes
  a conservative approach and complains that the featureCounts table does not
  match the sampletable. Specify ``subset.counts=TRUE`` to indicate that it's
  OK.
- the ``combine.by`` is applied to everything in the list; in this example, all
  counts for lines in the sample table that share the same "biorep" value will
  be summed.
- the ``parallel`` argument is also used for everything in the list

See the help for ``lcdbwf::make.dds`` for more details.

This chunk becomes a dependency of all of the ``results`` chunks below.

``results``
-----------

This is actually a series of chunks where the bulk of the differential
expression analysis takes place.

For simple cases, you probably just need one of these. But for complex
experimental designs where you end up doing lots of contrasts, it can get time
consuming to run them every time you change the RMarkdown file.

The end result of these chunks is a single list containing DESeq2 results
objects and metadata in (sub)lists. Each of these sublists has:

- ``res``, the results object
- ``dds``, the string name in ``names(dds.list)``
- ``label``, a "nice" label to use

A two-contrast list might look like this. This continues our example from above,
where we want to run the same contrast on all samples and after removing
replicate 4:

.. code-block:: r

    res.list <- list(

        # First contrast using all samples
        ko.vs.wt=list(
            res=results(
                dds.list[["main"]],
                contrast=c("genotype", "KO", "WT"),
                parallel=parallel
            ),
            dds=dds.list[["main""]],
            label="KO vs WT"
        ),

        # Same contrast, but use the dds object that had replicate 4 removed
        ko.vs.wt.no.rep4=list(
            res=results(
                dds.list[["no.rep.4"]],
                contrast=c("genotype", "KO", "WT"),
                parallel=parallel
            ),
            dds=dds.list[["no.rep.4""]],
            label="KO vs WT, without replicate 4"
        )
    )

If you have a small number of contrasts, this works fine. For complex
experimental designs, read on....

Complex experimental designs with many contrasts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For complex experimental designs with many contrasts, we can take advantage of
the ``knitr`` package's caching functionality to incrementally build results
objects and cache them. However, if we put all ``results()`` calls into the same
chunk and cache that, then a change anywhere in that chunk will invalidate the
cache, causing it all to be run again. An alternative is to put each
``results()`` call into its own chunk. But then we need to keep track of
dependencies and ensure those dependencies are specified in downstream chunks.
If you add a chunk and cache it, but forget to add the dependency later, the
environment will be inconsistent.

We use a modification of this second strategy. In the example above where we
have two contrasts (labeled ``ko.vs.wt`` and ``ko.vs.wt.no.rep4``), we put each
of those goes into its own chunk, but according to the following rules:

- the chunk name must start with ``results_``
- the variable name starts with ``contr_``, and the rest of the variable name
  will be used as the name in the list

So the above example becomes:

.. code-block:: r

    ```{r results_01, cache=TRUE, dependson=c('dds.list')}
        # First contrast using all samples
        contr_ko.vs.wt <- list(
            res=results(
                dds.list[["main"]],
                contrast=c("genotype", "KO", "WT"),
                parallel=parallel
            ),
            dds=dds.list[["main""]],
            label="KO vs WT"
        )
    ```

    ```{r results_02, cache=TRUE, dependson=c('dds.list')}
        # Same contrast, but use the dds object that had replicate 4 removed
        contr_ko.vs.wt.no.rep4 <- list(
            res=results(
                dds.list[["no.rep.4"]],
                contrast=c("genotype", "KO", "WT"),
                parallel=parallel
            ),
            dds=dds.list[["no.rep.4""]],
            label="KO vs WT, without replicate 4"
        )
    ```

Then we assemble everything together in a later chunk. The first trick in this
assembly chunk is that, because of the chunk naming scheme (names starting with
``results_``), we can automatically compile the list of chunks that are
dependencies. This ensures that the assembled list is up-to-date. The second
trick is that it inspects the environment to find variables with the naming
scheme ``contr_`` and does the work of inserting them into a list where the
names of the list come from the variable names without the ``contr_`` prefix.

.. code-block:: r

   ```{r assemble_variables, cache=TRUE, dependson=knitr::all_labels()[grepl('^results', knitr::all_labels())]}
    res.list <- list()
    contrast_list <- ls()[grepl("^contr_", ls())]
    res.list <- map(contrast_list, function(x) eval(parse(text=x)))
    res_names <- map(contrast_list, function(x) str_replace(x, "contr_", ""))
    names(res.list) <- res_names

The end result of this chunk is a list of lists that is used by functions in
the `lcdbwf` R package for more downstream work. For more details, see
:term:`res.list`.

For each contrast (that is, each entry in `res.list`) the below chunks will
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


Notes on using lfcShrink
^^^^^^^^^^^^^^^^^^^^^^^^
As currently implemented (05 apr 2018), lfcShrink checks its arguments for an
existing results table. If it exists, it applies shrinkage to the lfc and se
in that table. If it *doesn't* exist, it calls results on dds with the syntax

    res <- results(dds, name=coef)

or

    res <- results(dds, contrast=contrast)

It does not pass any further arguments to results, and it doesn't warn you
that results-style arguments were unrecognized and ignored. Therefore,
lfcShrink DOES NOT directly support lfcThreshold, or other alternative
hypotheses, or any of the custom analysis methods you can access through
results(). To get those, you have to call results first, without shrinkage,
and then apply lfcShrink.

Here we use the lfcShrink version of the results. In DESeq2 versions >1.16,
the lfc shrinkage is performed in a separate step, so that's what we do here.
This is slightly different results than if you used betaPrior=TRUE when
creating the DESeq object.



``attach``
----------

Typically the genes as labeled in the counts tables use Ensembl or other
not-quite-human-readable names. This chunk allows you to add additional gene
information to the results objects.

+---------+----------------------------------------------------------------------------------+
| var     | description                                                                      |
+=========+==================================================================================+
| keytype | in the counts table, what format are the gene IDs? Must be a column in the OrgDb |
+---------+----------------------------------------------------------------------------------+
| columns | what additional gene IDs to add? Must be columns in the OrgDb                    |
+---------+----------------------------------------------------------------------------------+

``reportresults``
-----------------

This is the section that creates multiple, tabbed outputs for each of the
contrasts in the :term:`res.list`.

``selections``
--------------

Here we get a list of DE genes from the :term:`res.list` object
to use for downstream analysis using the log2FoldChange (lfc) and 
false discovery rate (FDR) thresholds.

``upsetplots``
--------------

This chunk produces Upset plots comparing the selected lists of genes.

``helpdocs``
------------

For new users, or when distributing the output to collaborators who might
not be familiar with the plots contained in the report, a background and
help section are included as a child Rmd. This can be disabled by setting
`eval=FALSE` for this chunk.

``genepatterns``
----------------

Here we perform pattern analysis of the differentially expressed genes
to find co-regulated sets of genes using the ``DEGreport`` R package
in a separate child Rmd. For more details see `Gene patterns analysis`_ below.

``functionalenrichment``
------------------------

Here we perform functional enrichment analysis of the differentially expressed genes
to find enriched functional terms or pathways using the ``clusterProfiler`` R package.
This analysis is also performed in a separate child Rmd; for more details see `Functional enrichment analysis`_ below.

Gene patterns analysis
~~~~~~~~~~~~~~~~~~~~~~
See :ref:`gene-patterns` for details.

Functional enrichment analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
See :ref:`functional-enrichment` for details.

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

   res.list
      A list, with one item per contrast. Each of those items in turn is a list
      of objects that together compose the contrast (dds, results object, and
      label). This list-of-lists, which we call `res.list` for short, is used
      by functions in the `lcdbwf` R package for more downstream work.

      For a single contrast, it might look something like this:

      .. code-block:: r

         res.list[['contrast1']][['dds']] <- dds
         res.list[['contrast1']][['res']] <- res
         res.list[['contrast1']][['label']] <- 'Treatment vs control'

