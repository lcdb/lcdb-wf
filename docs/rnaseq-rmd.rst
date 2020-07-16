``rnaseq.Rmd`` code documentation
=================================


This page serves as documentation for the
:file:`workflows/rnaseq/downstream/rnaseq.Rmd` file, allowing the full power of
reStructured Text for presentation.

It is organized by chunk. For example, if :file:`rnaseq.Rmd` has the following
code:

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
documented and all documentation corresponds to a chunk.


``global_options``
------------------
This chunk sets global rmarkdown options. The following lines provide
a mechanism for all cached chunks to have their cache invalidated if either of
the filenames' modification times have changed.

.. code-block:: r

    cache.extra_file_dep_1=file.info('../config/sampletable.tsv')$mtime,
    cache.extra_file_dep_2 = file.info('../data/rnaseq_aggregation/featurecounts.txt')$mtime


``load_helpers``
----------------

The helper functions are stored in a package that is maintained in this
repository, in the ``lib/lcdbwf`` R package.

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

This is a good place to put any modifications to the sample table (like factors
derived other columns).

If you are using this Rmd outside the context of lcdb-wf, you will need to
change the salmon output path patterns and the sampletable location.

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



Note on factors
~~~~~~~~~~~~~~~
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
settled on the approach of a named list of ``dds`` objects.

**The** ``results`` **chunk below expects a list, one item per** ``dds`` **object.**

The simplest example is the following where we create a single ``dds`` and put
it into a list.

.. code-block:: r

   dds <- DESeqFromCombinedFeatureCounts(
      '../data/rnaseq_aggregation/featurecounts.txt',
      sampletable=colData,
      design=~group)
   dds <- DESeq(dds, parallel=parallel)

   dds.list <- list(main=dds)


Here is a modified example where we now want to remove replicate 4. We also want to collapse technical replicates:

.. code-block:: r

   dds1 <- DESeqFromCombinedFeatureCounts(
      '../data/rnaseq_aggregation/featurecounts.txt',
      sampletable=colData,
      design=~group)
   dds1 <- collapseReplicates(dds1, 'biorep')
   dds1 <- DESeq(dds1, parallel=parallel)

   dds2 <- DESeqFromCombinedFeatureCounts(
      '../data/rnaseq_aggregation/featurecounts.txt',
      sampletable=colData %>% filter(replicate!='rep4'),
      design=~group,
      subset.counts=TRUE  # need this to subset the featureCounts to match the colData
      )
   dds2 <- collapseReplicates(dds, 'biorep')
   dds2 <- DESeq(dds2, parallel=parallel)

   dds.list <- list(main=dds1, no.rep.4=dds2)

Based on our experience, as we add more ``dds`` objects the code gets more
error-prone. So for more complex use-cases, we have a function
``lcdbwf::make.dds``. This takes as its first argument a list of sampletable
(:term:`colData`) and a design and additional arguments can configure the
object further.

The above example becomes the following:

.. code-block:: r

   lst <- list(
      main=list(sampletable=colData, design=~group),
      no.rep.4=list(
         sampletable=colData %>% filter(replicate!='rep4),
         design=~group,
         args=list(subset.counts=TRUE))
   )

   dds.list <- map(lst, make.dds, combine.by='biorep', parallel=parallel)

Note the following:

- the file is set by default to be :file:`../data/rnaseq_aggregation/featurecounts.txt`
- we can supply additional args, like ``subset.counts=TRUE``, on a per-``dds`` basis.
- the `combine.by` is applied to everything in the list
- the ``parallel`` argument is also used for everything in the list

See the help for ``lcdbwf::make.dds`` for more details.

``results``
-----------

This chunk is where the bulk of the differential expression analysis takes place.

The end result of this chunk is a list of listes that is used by functions in
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


NOTE: Here are some notes on using lfcShrink...
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

res.list is a named list. Each item should be a list with names c('res',
'dds', 'label'). "res" is a DESeqResults object, "dds" is the corresponding
DESeq object the results were extracted from, and "label" is a nicer label to
use for headers and other text.


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



``upsetplots``
--------------

``helpdocs``
------------

``genepatterns``
----------------

``functionalenrichment``
------------------------


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

