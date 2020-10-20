.. _functional-enrichment:

Functional enrichment analysis
==============================

Let's say that the differential expression analysis of your RNA-Seq data
yields hundreds of changed genes. While this is exciting as it shows that
there is an effect of the treatment or mutant or drug, it also can be
daunting, since it is extremely labor intensive to manually go through a
bunch of genes, and assess their relevance to the biology that you're 
interested in. 

So, one of the first questions we usually have after performing a
differential expression analysis, is: **what kind of genes are changing in 
the experiment?**

To answer this question, we perform functional enrichment analysis, where
information about gene function (``annotation``) from various databases
are collated and tested for enrichment using a `hypergeometric test 
<http://en.wikipedia.org/wiki/Hypergeometric_distribution#Hypergeometric_test>`_.
The databases implemented here by default are gene ontology (GO) and KEGG pathways,
while Reactome pathways can be enabled if desired.

cprofsetup
----------
Here we set up ``clusterProfiler`` parameters and output directories. The default
parameters should work with almost all analyses. The only parameter that can be
changed:

+-----------+--------------------------------------------------------------------------------------------------+
| parameter | description                                                                                      |
+===========+==================================================================================================+
| ``types`` | A list of keys to extract from the annotation database. The defaults should not be changed, but  |
|           | additional keys can be added. For a list of keys for an orgDb object, run ``keyTypes(orgdb)``    |
+-----------+--------------------------------------------------------------------------------------------------+

``bitrgenes``
-------------
Here we obtain a list of differentially expressed genes and obtain a mapping of
these gene IDs to the terms listed in ``types`` above. Notable parameter to change:

+---------------+-----------------------------------------------------------------------+
| parameter     | description                                                           |
+===============+=======================================================================+
| ``from.type`` | This needs to be changed to match the source of gene IDs being used.  |
|               | For instance, if the GTF file is downloaded from Ensembl, then this   |
|               | should be set to ``ENSEMBL`` (default = ``'FLYBASE'``).               |
+---------------+-----------------------------------------------------------------------+

``enrichall``
-------------
This is the main chunk where the functional enrichment analysis is performed. A set of 
parameters for the analysis are set at the top of the chunk, but the only ones that should be
changed based on the experiment are:

+-----------------------+-----------------------------------------------------------------------------------+
| parameter             | description                                                                       |
+=======================+===================================================================================+
| ``kegg.organism``     | This is set based on the organism being studied. For example, for *Homo sapiens*  |
|                       | this should be ``'hsa'`` (default = ``'dme'``).                                   |
+-----------------------+-----------------------------------------------------------------------------------+
| ``RUN.REACTOME``      | Specify if Reactome analysis should be performed (default = ``FALSE``).           |
+-----------------------+-----------------------------------------------------------------------------------+
| ``reactome.organism`` | In contrast to KEGG, Reactome analysis needs the name of the species.             |
|                       | For example, for *Homo sapiens* this should be ``'human'`` (default = ``'fly'``). |
|                       | This only needs adjusting if ``RUN.REACTOME`` is ``TRUE``.                        |
+-----------------------+-----------------------------------------------------------------------------------+

This chunk performs GO, KEGG (and optionally Reactome) analyses for each element in the list of 
DE genes, separately for up- and down-regulated genes. The results populate a list which is then
output to files and used subsequently for visualization.

.. warning::

   Each list of DE genes is analyzed in multiple (3 GO ontologies + KEGG + optionally Reactome) separate ways
   so for many lists, this chunk can be *very* time-consuming. However, with ``cache = TRUE``, 
   this is only run once and subsequent runs are much faster.

``plotgo``
----------
Here we visualize the functional enrichment analysis results using three different
plots:

- ``dotplot``

   This represents highly enriched terms from the analysis with the color corresponding to
   level of enrichment (adjusted p-value) and the size of the dots representing the
   number of genes associated with the GO term or pathway.
   
   This is the easiest plot to read the enriched terms. By default the top 10
   categories are shown.

- ``emapplot``

   The enrichment map plot is a network visualization of the enrichment analysis. Nodes in
   the network represent GO terms or pathways and the connections (edges) between nodes
   represent shared genes between the terms/pathways. The color of the node corresponds
   to adjusted p-values, the size represent the number of genes associated with the term
   and the thickness of the edge represents the number or fraction of genes shared between
   the terms.
   
   This plot is good for identifying enriched pathways, since categories with
   similar genes will cluster together. By default the top 30 categories are
   shown.

- ``cnetplot``

   This plot extracts the complex associations of genes with multiple functional terms or
   pathways. The size of the nodes for functional terms represent the number of associated
   genes and the gene nodes can also be overlayed with observed fold-changes.
   
   This plot is good for identifying which genes act as links between which
   categories. By default the top 5 categories are plotted.

Each plot is also saved to file as PDF and a link included in the html summary report.

``adjustcatlen``
----------------
Here we update the ``clusterProfiler`` results object by adding symbols corresponding to
KEGG gene IDs. This is done so that human-readable gene names show in the visualizations.
This will likely be moved to a helper function that is run automatically in a future release.

