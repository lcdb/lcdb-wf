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


global_options
--------------
The following lines provide a mechanism for all cached chunks to have their
cache invalidated if either of the filenames' modification times have changed.

.. code-block:: r

    cache.extra_file_dep_1=file.info('../config/sampletable.tsv')$mtime,
    cache.extra_file_dep_2 = file.info('../data/rnaseq_aggregation/featurecounts.txt')$mtime

imports
-------


run_helpers
-----------




annotationhub_setup
-------------------

We use AnnotationHub for downloading annotations (rather than installing
various OrgDb packages). Once downloaded, this is cached in the
:file:`../../../include/AnnotationHubCache` directory so it can be used by
other workflows.

Specify the g

coldata_setup
-------------

salmon
------

ddstxi
------

dds_initial
-----------

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
