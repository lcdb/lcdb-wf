lcdb-wf
=======

`lcdb-wf` is a collection of Snakemake workflows for common high-throughput
sequencing analysis.

There are a lot of workflows out there for high-throughput sequencing analysis.
What makes `lcdb-wf` different?

Designed with customization in mind
-----------------------------------
We recognize that every experiment has its own unique idiosyncracies. Rather
than provide a one-size-fits-all solution, we aim to provide a reasonable
starting point that users can modify for their own use.

Non-model organism? Custom gene annotations? Complicated regression models?
Unconventional command-line arguments to tools? New tools to add to the
workflow? No problem.

Tested automatically
--------------------
Every change to the code on GitHub triggers an automated test, the results of
which you can find at https://circleci.com/gh/lcdb/lcdb-wf. Each test sets the
system up from scratch, including installing all software, downloading example
data, and running everything up through the final results.

Track hubs
----------
The ChIP-seq and RNA-seq workflows generate track hubs that can be viewed in
the UCSC Genome Browser. ChIP-seq shows signal and called peaks (across
multiple peak-callers); RNA-seq shows strand-specific signal tracks. Both
support the addition of arbitrary additional tracks (primers, loci of interest,
external data, etc) to view alongside your data.

Unified approach to references
------------------------------
The references workflow is shared by RNA-seq and ChIP-seq and is driven by
a config file that specifies URLs for FASTA and GTF files. Set it up once for
a site to get lots of genomes you can use for running `fastq_screen`, and
easily include arbitrary other genomes. They can then be automatically included
in RNA-seq and ChIP-seq workflows.

Integration with external data and figure-making
------------------------------------------------
Examples are included (and tested) for how to tie together workflows, including
downloading external data and downstream figure-making work.

If an upstream file changes (e.g., gene annotation), all downstream jobs
depending on it -- including figures -- will be updated so you can ensure that
even complex analyses stay correct and up-to-date.

Advantages of Snakemake
-----------------------

The workflows are built using `Snakemake
<https://snakemake.readthedocs.io/en/stable/>`_, so we also get the following
for free:

Run locally or on a cluster with the same code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Run the same workflow locally or on a cluster with only a slight command-line
modification. SLURM support for the NIH Biowulf cluster is built in; edit the
``clusterconfig.yaml`` file to specify cluster arguments specific to your
system.

Only run the jobs required
~~~~~~~~~~~~~~~~~~~~~~~~~~
New gene annotation? Snakemake tracks dependencies, so it will detect that the
annotations changed, and re-run only the jobs that depend on that file at some
point in their dependency chain, leaving other files that don't depend on it
untouched.

Parallelization
~~~~~~~~~~~~~~~
It's trivial to run jobs in parallel, for as many CPUs or cluster nodes as you
have available by using the ``-j`` argument to Snakemake.

Software installation
~~~~~~~~~~~~~~~~~~~~~
Installation of all dependencies is handled by conda, ensuring reproducibility,
streamlined setup, and no need for root administrator privileges.

See :ref:`getting-started` to get started.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   getting-started
   tests
   workflows
   config
   references
   rnaseq
   downstream-rnaseq
   chipseq
   external
   figures
   colocalization
   wrappers
   cluster
   troubleshooting
   autodoc
   changelog
   developers
