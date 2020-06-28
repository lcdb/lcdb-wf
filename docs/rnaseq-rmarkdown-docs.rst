``rnaseq.Rmd`` code documentation
=================================

Previous versions included extensive documentation within the code, but these
comments got too large and were still incomplete. This page serves as the
documentation instead, allowing the full power of reStructured Text for
presentation.

It is organized by chunk. For example, if the rnaseq.Rmd has the following code:

.. code-block:: text

    ```{r load_libraries}
    library(DESeq2)
    library(dplyr)
    ```

Then you can find the corresponding documentation on this page under the
``load_libraries`` heading. Note that an RMarkdown code chunk must be named in
order to be documented. RMarkdown requires uniquely-named chunks, so we can use
them as uniquely-named headings.


.. note::

    Some terminology:

    :dds:
        DESeq data set

    :rld:
        The variance-stabilized transformed version of the counts. Used for
        PCA, clustered heatmaps, and gene patterns.

global_options
--------------
The following lines provide a mechanism for all cached chunks to have their
cache invalidated if either of the filenames' modification times have changed.

.. code-block:: r

    cache.extra_file_dep_1=file.info('../config/sampletable.tsv')$mtime,
    cache.extra_file_dep_2 = file.info('../data/rnaseq_aggregation/featurecounts.txt')$mtime

imports
-------

All necessary libraries are loaded here.


run_helpers
-----------

This is the first of several [child
chunks](https://bookdown.org/yihui/rmarkdown-cookbook/child-document.html). The
main content lives in another Rmd and it is included directly into this Rmd.

When rendering the entire document at once, the contents of the chunk are
skipped. However, when interactively running chunk-by-chunk, there is no good
mechanism for this. Therefore, a line is included in the chunk to run the child
Rmd -- but don't run pandoc on it to get a separate HTML. The end result is
that the child is always run whether you're rendering the whole document or
just running chunk-by-chunk.


annotationhub_setup
-------------------

We use AnnotationHub for downloading annotations (rather than installing
various OrgDb packages). Once downloaded, this is cached in the
:file:`../../../include/AnnotationHubCache` directory so it can be used by
other workflows.

Providing the species to the `get_orgdb` helper function will search for the
species and return the latest OrgDb available in AnnotationHub.

coldata_setup
-------------

This chunk prepares the colData object that contains the metadata to be used
for creating the dds objects.

The majority of the metadata is read in from the sampletable. This chunk also
handles things like converting to factors and setting up the paths to Salmon
output files.

This is a good place to put any modifications to the sample table (like factors
derived other columns).

If you are using this Rmd outside the context of lcdb-wf, you will need to
change the salmon output path patterns and the sampletable location.

salmon
------

If you don't want to use Salmon TPM, disable this chunk with ``eval=FALSE`` or
delete it entirely (and do the same with the next chunk).

ddstxi
------

This creates separate ``dds.txi`` and ``rld.txi`` objects to differentiate them
from the ones with no ``.txi`` that are created using featureCounts.

dds_initial
-----------
Create the initial dds and rld objects that will be used for exploratory data analysis.


dds_models
----------

results
-------

attach
------

selections
----------

upsetplots
----------

helpdocs
--------

child='gene-patterns.Rmd'
-------------------------

child='functional-enrichment.Rmd'
---------------------------------
