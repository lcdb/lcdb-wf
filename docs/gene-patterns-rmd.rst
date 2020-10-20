.. _gene-patterns:

Gene patterns analysis
======================

When analyzing RNA-Seq experiments that have a temporal element, e.g.
time-series data where a treatment is applied for an increasing amount
of time, or a dose-response aspect, e.g. where increasing amounts of a drug
is administered to a set of cells, an interesting biological question can be:
*which genes are changing in the same way over time, or in response to the
drug treatment?*

In these cases, it may be useful to look for groups of **co-expressed** genes, or genes
that show similar trends in expression over time or the dose-response assay, as this might
indicate a pathway or gene regulatory network is being affected. For instance,
if the main effect of a treatment is to suppress a transcription factor that
induces a group of genes, these genes would likely show a trend of decreasing
expression over the course of a dose-response experiment.

One approach to this analysis could be to compare each time-point or dose sample to
the control sample. However, if we are working with many time-points or doses
this can get unwieldy very quickly, as we would be dealing with that many lists
of differentially expressed genes. Moreover, we would lose the temporal information.
For instance, in a time-course experiment if a gene is up-regulated from 0 hr ⇒ 2 hr, 
and also from 0 hr ⇒ 4 hr, we don't know how 4 hr compares to 2 hr, without comparing
those time-points explicitly.

An alternative approach that avoids the complications of a per-time-point analysis
is a pattern-based analysis, which looks at the temporal trend of gene expression as a whole,
to find similar groups of genes. For instance, in the above experiment, the pattern
analysis would immediately tell us which genes keep going up from 0 ⇒ 2 ⇒ 4 hr,
and which genes go up from 0 ⇒ 2 hr, but then don't change from 2 ⇒ 4 hr.

In this Rmd, we implement the latter approach and look for groups of co-expressed or
co-regulated genes, using gene patterns analysis with the R package, ``DEGreport``. The input to 
this analysis is normalized expression data of differentially expressed genes.
A clustering algorithm implemented in the package is then used to find similar groups
of genes, which are then plotted.

``dpsettings``
--------------

This is the chunk where the parameters for the pattern analysis are specified.
Notable parameters to adjust include:

+-------------+----------------------------------------------------------------------------------------------------------------+
| parameter   |  description                                                                                                   |
+=============+================================================================================================================+
| ``time``    | Factor to show on x-axis. Typically this is the time or dose-response factor in the experiment.                |
|             | Should be a column in the ``colData`` metadata (default = ``'group'``).                                        | 
+-------------+----------------------------------------------------------------------------------------------------------------+
| ``col``     | Factor to color lines by. Should be a column in the ``colData`` metadata (default = ``NULL``).                 |
+-------------+----------------------------------------------------------------------------------------------------------------+
| ``minc``    | Minimum cluster size (default = 1). Consider increasing if you're getting many clusters with very few genes.   |
+-------------+----------------------------------------------------------------------------------------------------------------+
| ``lim``     | Maximum number of genes to include in the clustering (default = 2000). If number of DE genes is higher, they   |
|             | are downsampled to ``lim``. Increase with caution since clustering a huge number of genes can be very          | 
|             | CPU-intensive.                                                                                                 |
+-------------+----------------------------------------------------------------------------------------------------------------+

``finalclusters``
-----------------

This is the chunk where the pattern analysis is done via the following steps:

1. Set up selections of DE genes to be used in the clustering analysis. By default,
   we use all the ``changed`` genes for this purpose, but you can choose to only
   use ``up`` or ``dn`` genes.
2. The set of genes are consolidated and downsampled if the total number > ``lim``.
3. Normalized data corresponding to these genes are extracted. We use the
   ``varianceStabilizingTransformation`` for the normalization. Before proceeding,
   genes with identical counts across all samples are removed, since this could lead
   to errors in the clustering.
4. The primary clustering command ``degPatterns`` is run, using the parameters specified
   in the ``dpsettings`` chunk above.
5. Clusters having > ``minc`` genes are filtered out and plotted, with the gene lists being
   saved to files in the ``final_clusters`` subdirectory, while a column encoding the cluster
   membership is added to the ``res.list`` element corresponding to the DE genes being analyzed.
