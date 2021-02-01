FAQs
====

This page serves as a catch-all for details of various topics.


.. _simultaneous-workflows:

Can I run multiple workflows at once?
-------------------------------------

Sometimes. While Snakemake creates a lockfile to prevent multiple instances
running in the same directory, the fact that each workflow (rnaseq, chipseq,
etc) are in their own subdirectory means that they will each have their own
separate lockfiles and can be run.

Be careful though, because both ChIP-seq and RNA-seq workflows include the
references workflow. This means that if you have not yet already set up the
references, the RNA-seq and ChIP-seq workflows may both attempt to write the
references, potentially corrupting it.


.. _multiple-experiments:

How do I handle multiple experiments in the same project?
---------------------------------------------------------

It's pretty common to have both RNA-seq and ChIP-seq experiments that need to
be analyzed together. For example, we might have RNA-seq in two different mouse
cell types, RNA-seq in a human cell type, ChIP-seq for different antibodies in
those cell types, and all of this needs to be compared with publicly available
data (say, other RNA-seq experiments from GEO). We need to make figures for the
manuscript, so the figure-making code needs to be re-run if there are any
changes to the primary analysis (new samples, parameter changes, etc).

lcdb-wf is designed to handle all of this. There are a couple of limitations to
lcdb-wf that will determine how best to split up workflows, and the subsections
below will help you decide if you should consider an experiment as part of
a different workflow.

If an experiment needs to be considered as part of a different workflow, then
make copies of the relevant workflow directory after deploying. Taking the
above project as an example, immediately after deploying (using ``--flavor
full`` so we get all supported workflows including RNA-seq and ChIP-seq), we
have::

    workflows/
      chipseq/
      colocalization/
      external/
      figures/
      rnaseq/

then we might rename the directory called "external" to match the GEO accession (to make it
easier to remember), make copies of the RNA-seq directory for the mouse and
human experiments, and clean up a little:

.. code-block:: bash

    rm -r workflows/colocalization
    mv workflows/external workflows/GSE00112233
    cp -r workflows/rnaseq workflows/mouse-rnaseq
    cp -r workflows/rnaseq workflows/human-rnaseq
    rm -r workflows/rnaseq

then we would have::

    workflows/
      chipseq/
      figures/
      GSE00112233/
      human-rnaseq/
      mouse-rnaseq/

See below for advice on when to split experiments into separate workflows.


Samples need the same library layout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A single workflow must be either all single-end or all paired-end.

If you have the same RNA-seq samples that were sequenced once SE and again PE,
you'll need to have two different copies of the ``workflows/rnaseq`` directory
(say, ``workflows/rnaseq-se`` and ``workflows/rnaseq-pe``). If they are to be
combined in the same differential expression analysis (e.g., by modeling layout
as a batch effect), then the adjust your downstream analysis in R to read in
both counts tables (see :ref:`rnaseqrmd` for more).

Samples need to use the same parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is no mechanism for specifying sample-specific parameters. For example,
to use cutadapt to trim some 5' bases from some samples. Other samples are left
alone. Samples that need to treated differently should be split off into
a separate workflow, and the respective Snakefiles should be edited
accordingly.

.. note::

    Note that the peak-calling for ChIP-seq supports specifying custom
    parameters for each peak-calling run. The BAM files still need to have used
    uniform parameters across samples, so if alignment or trimming options
    differ, you should split each set of parameters into a different workflow.

Samples must use the same assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is no mechanism for specifying per-sample assemblies. Samples from
different species need to go in different workflows. If you want to compare,
say, hg19 and hg38, then you would make a copy of the workflow for each assembly
and adjust the reference config accordingly for each workflow separately.

.. _lowcounts:

How are low counts handled during differential expression analysis? Should we use a read-count threshold to filter genes?
-------------------------------------------------------------------------------------------------------------------------
Low count genes are handled during the normalization and analysis steps of DESeq2
with sophisticated statistical models. Genes with low counts across the board are flagged
as *low count outliers*, and the p-values are set to NA. Also genes with low counts
are penalized by shrinking the ``log2FoldChange`` estimate. For example, a fold change of
4 that came from 4 reads in the treatment group vs 1 read in the control, will be shrunken,
as opposed to if the treatment had 2000 reads vs 500 in the control. As a result of this
low-count correction, the ``log2FoldChange`` of genes clearing a false-discovery criterion
can be used as a reliable metric for prioritizing candidate genes for follow-up experiments.
In contrast, using an arbitrary fold-change cutoff could introduce biases that potentially
violate modeling assumptions and introduce variables that we could not predict or control for.
So, we do not recommend using count thresholds to filter differential expression analysis
results to determine candidate genes for follow up.


.. _troubleshooting:

How do I troubleshoot failed jobs?
----------------------------------
Many rules have an explicit ``log:`` directive that defines where the log is
written. These are typically in the same directory as the output files the rule
creates, and this is the first place to check if something goes wrong.

Some rules do not explicitly redirect to ``log:`` or may only redirect either
stdout or stderr. Where this output ends up depends on if you're running
locally or on a cluster.

**When running locally,**  stdout and stderr will be included in the output
from Snakemake, so check there.

**If running on a cluster,** the default behavior is to send the main Snakemake
output to ``Snakefile.log``.  The per-rule output depends on how it was sent to
the cluster.  As described in the above section, by default stdout and stderr
are sent to the ``logs`` directory, named after rule and job ID.

**If a job fails on a cluster**:

- Open ``Snakefile.log`` and search for ``Error``
- Recent versions of Snakemake report the ``log:`` file (if any) and the
  ``cluster_jobid:``. Keep track of these.
- If ``log:`` was defined for the rule, check there first
- If not, or if more information is needed, check
  ``logs/<rulename>.{e,o}.<jobid>`` (which is how stderr and stdout are
  configure when running with the ``include/WRAPPER_SLURM`` wrapper).

For example, if we find the following error in ``Snakefile.log``::

    [Tue Feb  6 20:06:30 2018] Error in rule rnaseq_rmarkdown:
    [Tue Feb  6 20:06:30 2018]     jobid: 156
    [Tue Feb  6 20:06:30 2018]     output: downstream/rnaseq.html
    [Tue Feb  6 20:06:30 2018]     cluster_jobid: 60894387

Then we would check ``logs/rnaseq_markdown.e.60894387`` and
``logs/rnaseq_markdown.o.60894387`` for more information.


.. _updating:

How do I update my deployment?
------------------------------

If there are additional fixes or features in the main lcdb-wf repo that you
want to propagate to your existing projects, the best way to do this is to
clone a recent version and do the manual diffs between the new version and what
you have on disk.

To help narrow down the changes that have happened in the main lcdb-wf repo
since you deplyed to a project, Use the ``.lcdb-wf-deployment.json`` file that
is created when deploying to a project to find the commit hash that the
deployment used.
