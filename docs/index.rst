Introduction
============

`lcdb-wf` is a collection of Snakemake workflows for common high-throughput
sequencing analysis.

There are a multitude of workflows out there for high-throughput sequencing analysis.
What makes `lcdb-wf` different?

Designed with customization in mind
-----------------------------------
We recognize that every experiment has its own idiosyncracies. Rather
than provide a one-size-fits-all solution, we aim to provide a reasonable
starting point that users can modify for their own use.

Non-model organism? Custom gene annotations? Complicated regression models?
Unconventional command-line arguments to tools? New tools to add to the
workflow? No problem.

Extensive downstream RNA-seq
----------------------------
A comprehensive RMarkdown template, along with a custom R package, enables
sophisticated RNA-seq analysis that supports complex experimental designs and
many contrasts.

Extenstive exploration of ChIP-seq peaks
----------------------------------------
The ChIP-seq configuration supports multiple peak-callers as well as calling
peaks with many different parameter sets for each caller. Combined with
visualizaiton in track hubs (see below), this can identify the optimal
parameters for a given experiment.

Track hubs
----------
The ChIP-seq and RNA-seq workflows generate track hubs that can be viewed in
the UCSC Genome Browser. ChIP-seq shows signal and called peaks (many
peak-calling runs, with different peak-callers and different parameters for
each, can be configured to automatically run and show up in the track hub);
RNA-seq shows strand-specific signal tracks. Both support the addition of
arbitrary additional tracks (primers, loci of interest, external data, etc) to
view alongside your data.

Unified approach to reference genomes
-------------------------------------
The references workflow defines the genome files for the organisms used in 
the experiment. It it is shared by RNA-seq and ChIP-seq and is driven by
a config file that specifies URLs for FASTA and GTF files. Set it up once for
a site to get lots of genomes you can use for running `fastq_screen`, and
easily include arbitrary other genomes. They can then be automatically included
in RNA-seq and ChIP-seq workflows.

Arbitrary genomes can be used, whether local (e.g., customized with additional
genetic constructs) or on the web. The `references` workflow need only be run
once for all these genomes to be created, with the `references_dir` being used
as a centralized repository that can be then used with all other workflows.

Integration with external data and figure-making
------------------------------------------------
It is designed to tie together your entire analysis, from downloading references
through producing final figures. The included examples tie together workflows that download external data, perform downstream analysis, and make figures.

If an upstream file changes (e.g., gene annotation), all dependent downstream
jobs -- including figures -- will be updated so you can ensure that even
complex analyses stay correct and up-to-date.

Tested automatically
--------------------
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
