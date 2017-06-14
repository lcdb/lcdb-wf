Features common to all workflows
================================
TMPDIR handling
---------------
The top of each snakefile sets up a shell prefix that exports the TMPDIR
variable. The reason for this is that the NIH Biowulf cluster supports nodes
with temporary local storage in a directory named after the SLURM job ID. This
ID is not known ahead of time, but is stored in the ``SLURM_JOBID`` env var.

Since each rule executed on a cluster node calls the snakefile (see the job
scripts created by snakemake for more on this), we can look for the job ID and
set the tempdir appropriately. Upon setting ``$TMPDIR``, the Python
``tempfile`` module will use that directory to store temp files. Any wrappers
can additionally use ``$TMPDIR`` in shell commands and it will use this
directory.

Note that the default behavior is to set ``$TMPDIR`` to the default temp
directory as documented in Python's `tempfile module
< https://docs.python.org/3/library/tempfile.html#tempfile.gettempdir>`_.

References dir
--------------
A references directory must be set either as an environment variable
``REFERENCES_DIR`` or in the config file in the ``references_dir`` key. This
specifies the top-level directory in which references will be built. See
:ref:`references` for more details.

Sample table and samples
------------------------
Most workflows (with the exception of ``references.snakefile``) use a sample
table. This is specified in the config file under the ``sampletable`` key. The
first column is assumed to be the sample IDs.

`patterns` and `targets`
------------------------
Workflows have a `patterns` dict at the top that lays out, in one place, the
output file patterns. Rules can optionally use values from this dict to define
output patterns. Using the `lcdblib.snakemake.fill_patterns` function, these
patterns are "rendered" into a `targets` dictionary (this is basically
a recursive `expand()` using rows from the sampletable). Selected contents of
the `targets` directory are then used for the `all` rule.

This provides, all in one section of the Snakefile, a fair amount of
customization. Patterns can be commented out, or targets can be excluded from
the `all` rule to fine-tune which rules will be run.

A nice side-effect of the `patterns` and `targets` design is that aggregation
rules become much easier to write. If all FastQC runs are under a `fastqc` key
in `targets`, they can all be used as input to a rule using
`lcdblib.utils.flatten(targets['fastqc'])`. This is in contrast to, say,
writing input functions to track down the various FastQC output files. Another
example of this is the inputs section of the ``multiqc`` rule in the RNA-seq
workflow.
