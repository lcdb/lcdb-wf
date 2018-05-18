.. _troubleshooting:

Troubleshooting
---------------
Many rules have an explicit ``log:`` directive that defines where the log is
written. These are typically in the same directory as the output files the rule
creates, and this is the first place to check if something goes wrong.

Some rules do not explicitly redirect to ``log:`` or may only redirect either
stdout or stderr. Where this output ends up depends on if you're running
locally or on a cluster. When running locally, stdout and stderr will be
included in the output from Snakemake, so check there.

If running on a cluster, the default behavior is to send the main Snakemake
output to ``Snakefile.log``.  The per-rule output depends on how it was sent to
the cluster.  As described in the above section, by default stdout and stderr
are sent to the ``logs`` directory, named after rule and job ID.

**If a job fails on a cluster**:

- Open ``Snakefile.log`` and search for ``Error``
- Recent versions of Snakemake report the ``log:`` file (if any) and the
  ``cluster_jobid:``. Keep track of these.
- If ``log:`` was defined, check there first for info
- If not, or if more information is needed, check
  ``logs/<rulename>.{e,o}.<jobid>``.

For example, if we find the following error in ``Snakefile.log``::

    [Tue Feb  6 20:06:30 2018] Error in rule rnaseq_rmarkdown:
    [Tue Feb  6 20:06:30 2018]     jobid: 156
    [Tue Feb  6 20:06:30 2018]     output: downstream/rnaseq.html
    [Tue Feb  6 20:06:30 2018]     cluster_jobid: 60894387

Then we would check ``logs/rnaseq_markdown.e.60894387`` and
``logs/rnaseq_markdown.o.60894387`` for more information.

