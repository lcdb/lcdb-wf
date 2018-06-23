
.. _colocalization:

Colocalization workflow
-----------------------
The output of this workflow is a set of heatmaps showing metrics of
colocalization between pairs of regions. This can be used to answer questions
like "what else does my protein of interest bind with?".

Several colocalization methods are run, and they all give slightly different
results.

- bedtools fisher
- bedtools jaccard
- GAT (Genome Association Test) log2 fold change
- GAT (Genome Association Test) nucleotide overlap
- IntervalStats
