Introduction
============

`lcdb-wf` is a collection of Snakemake workflows for common high-throughput
sequencing analysis.

What makes `lcdb-wf` different?

Features
--------
Designed with customization in mind
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We recognize that every experiment has its own idiosyncracies. Rather
than provide a one-size-fits-all solution, we aim to provide a reasonable
starting point that users can modify for their own use.

Non-model organism? Custom gene annotations? Complicated regression models?
Unconventional command-line arguments to tools? New tools to add to the
workflow? No problem.

Extensive downstream RNA-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A comprehensive RMarkdown template, along with a custom R package, enables
sophisticated RNA-seq analysis that supports complex experimental designs and
many contrasts. For example, easily set up multiple DESeqDataSet objects with
different models or different mixtures of samples, set up contrasts using these
different objects, and all the output and functional enrichment analysis will
automatically be generated.

Integration with Carnation
~~~~~~~~~~~~~~~~~~~~~~~~~~
The RNA-seq workflow generates RDS objects compatible with `Carnation
<https://github.com/NICHD-BSPC/carnation>`__ the Shiny app for interactively
exploring RNA-seq results, including gene patterns, functional enrichment,
comparisons across contrasts, and more.

Extensive exploration of ChIP-seq peaks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ChIP-seq configuration supports multiple peak-callers as well as calling
peaks with many different parameter sets for each caller. Combined with
visualizaiton in track hubs (see below), this can identify the optimal
parameters for a given experiment.

Track hubs
~~~~~~~~~~
The ChIP-seq and RNA-seq workflows generate track hubs that can be viewed in
the UCSC Genome Browser. ChIP-seq shows signal and called peaks (many
peak-calling runs, with different peak-callers and different parameters for
each, can be configured to automatically run and show up in the track hub);
RNA-seq shows strand-specific signal tracks. Both support the addition of
arbitrary additional tracks (primers, loci of interest, external data, etc) to
view alongside your data.

Support for complex reference genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Reference genomes may need to be patched with experimental genetic constructs,
or may need to be adjusted when downloaded from an original source (for
example, change chromosome nomenclature to match existing work).

Arbitrary genomes can be used, whether local (e.g., customized with additional
genetic constructs) or on the web.

Tested automatically
~~~~~~~~~~~~~~~~~~~~
Every change to the code on GitHub triggers an automated test, the results of
which you can find at https://circleci.com/gh/lcdb/lcdb-wf. Each test sets the
system up from scratch, including installing all software, downloading example
data, and running everything up through the final results. This guarantees that
you can set up and test the code yourself.

All the advantages of Snakemake
-------------------------------

The workflows are built using `Snakemake
<https://snakemake.readthedocs.io/en/stable/>`_, so we also get the following
for free:

Run locally or on a cluster with the same code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Run the same workflow locally or on a cluster with a single command-line flag.
Use `Snakemake profiles
<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ to
translate from general resources to cluster-specific
arguments.

Only run the required jobs
~~~~~~~~~~~~~~~~~~~~~~~~~~
New gene annotation? Snakemake tracks dependencies, so it will detect that the 
annotations changed. Only jobs that depend on that file at some point in their 
dependency chain will be re-run and the independent files are untouched. Adding
a new sample will leave untouched any output from samples that have already
run.

Parallelization
~~~~~~~~~~~~~~~
It's trivial to run jobs in parallel, for as many CPUs or cluster nodes as you
have available by using the ``-j`` argument to Snakemake.

Software installation
~~~~~~~~~~~~~~~~~~~~~
Installation of all dependencies is handled by conda, ensuring reproducibility,
streamlined setup, and no need for root administrator privileges.

See :ref:`getting-started` to get started.
