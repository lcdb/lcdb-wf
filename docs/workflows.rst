Features common to all workflows
================================

.. _config:


Configuration
-------------
In the top-level directory, you may want to edit ``include/WRAPPER_SLURM`` to
activate the environment you're using.

Within each workflow are the following files to check:

:``config/config.yaml``:
    This file contains information on what assembly to use, and configures the
    directories where output will go.

    Be sure to check:
    - assembly
    - references_dir
    - tags (for aligner, rrna, gtf, and salmon).

    You will probably want to copy-paste the references section from the
    top-level ``config/config.yaml``, since that has config information for
    many model organisms already.

:``Snakefile``:
    In the ``Snakefile``, search for the string ``TEST SETTINGS``.  This is used to
    indicate where settings can and should be customized, but are set by default to
    whatever makes sense for our test data. The comments often give alternative
    suggestions.

:``config/*_patterns.yaml``:
    The ChIP-seq and RNA-seq workflows have a "patterns" file. See below for
    more details on this, but generally unless you have strong opinions about
    output directory organization, you don't need to change this.

:``config/sampletable.tsv``:
    A tab-separated file containing metadata for the experiment. The first
    column is the sample name. The ``orig_filename`` is optional. By default,

:``config/clusterconfig.yaml``:
    Per-rule settings. Modify these settings according to your cluster. The
    existing settings are for NIH's Biowulf cluster, which runs the SLURM
    scheduler.

:``config/hub_config.yaml``:
    This file configures information about the track hub -- which colors to use
    for which samples, and user/hostname information for uploading the hub.


Running on a cluster
--------------------
The example commands in :ref:`getting-started` describe running Snakemake
locally. For larger data sets, you'll want to run them on an HPC cluster.
Snakemake `supports arbitrary cluster commands
<http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html>`_,
making it easy to run these workflows on many different cluster environments.

The default configuration we provide is specific to the NIH Biowulf cluster.
To run a workflow on Biowulf, from the workflow directory (e.g.,
``workflows/rnaseq``, run the following command::

    sbatch ../../include/WRAPPER_SLURM

The ``WRAPPER_SLURM`` script submits the main Snakemake process on a separate node to avoid any
restrictions from running on the head node. That main snakemake process then
submits each rule separately to the cluster scheduler. As configured in that
script, we specify ``config/clusterconfig.yaml`` as containing the
rule-specific cluster arguments.

That script also contains this Snakemake arguments::

    --jobname "s.{rulename}.{jobid}.sh" \
    --cluster 'sbatch {cluster.prefix} --cpus-per-task={threads}  --output=logs/{rule}.o.%j --error=logs/{rule}.e.%j' \

This means that each job is named after the rule and job id (the ``--jobname``
arg) and the stdout and stderr go to files in ``logs`` and are named after the
rule, followed by a ``.o`` or ``.e``, followed by the cluster job ID (the
``--cluster`` arg).


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


`patterns` and `targets`
------------------------
The input and output patterns are configured in a ``config/*_patterns.yaml``
file. This file provides a convenient way to organize the many files that are
created by a workflow. In general you can ignore this when using the workflows.
But if you want to create your own rules or otherwise do development work, it
may be helpful to understand this mechanism for the convenience it offers.

In each workflow, it is loaded into a ``SeqConfig``
object, which gives access to the patterns as well as the targets generated by
filling in those patterns with config information. For example in the RNA-seq
workflow, these patterns are in ``config/rnaseq_patterns.yaml``. They are
loaded into the RNASeqConfig object ``c``, and targets are filled in according
to the config (which includes all the sample names in the sampletable).
Patterns and targets for rules are accessed from this object.

To make it easier to access large lists of items, we can use the
``lcdblib.utils.flatten`` function to get a flat list of files starting from
any key in the patterns YAML file.

This has several advantages:

- Writing aggregation rules is much easier. For example, instead of lots of
  ``expand()`` calls, we can get all the FastQC output with
  ``flatten(c.targets["fastqc"])``.
- Toggling entire sections of the workflow can be performed by changing
  a single line in the first "all" rule. For example:

.. code-block:: python

    # include the dupradar files in the DAG
    rule all:
        input:
            flatten(c.targets['downstream']),
            flatten(c.targets['dupradar'])


    # exclude them
    rule all:
        input:
            flatten(c.targets['downstream']),
            # flatten(c.targets['dupradar'])

- You can re-organize the output directories by editing the patterns file
  without having to touch the Snakefile

- Downstream tasks can access these patterns and targets as well, avoiding the
  need to retype filename patterns and maintain them across multiple workflows.

To illustrate this last case, imagine we have the following analyses:

- RNA-seq data for a knockdown and control (in the ``rnaseq`` workflow)
- ChIP-seq data from some factors (analyzed with the ``chipseq`` workflow)
- some BED files of called peaks downloaded and lifted over from published data
  (``external`` workflow). 

Now we want to make a figure: an M-A plot of RNA-seq differential expression
with points colored by whether they have a ChIP-seq peak at their promoter. The
input for this figure comes from different workflows. It should be re-generated
if the differential expression changes, or if any of the peak calls change.

The way this is tyically handled is with Snakemake `subworkflows
<http://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#sub-workflows>`_,
and that's what we use here. However we can take advantage of the patterns and
targets to conveniently access any of the files created by any of the other
workflows. Here's an example:

.. code-block:: python

    # Create config objects. They need the config yaml, the patterns yaml, and
    # the working directory for the corresponding workflow.
    rnaseq_config = RNASeqConfig('config/config.yaml', 'config/rnaseq_patterns.yaml', workdir='../rnaseq')
    chipseq_config = ChIPSeqConfig('config/config.yaml', 'config/chipseq_patterns.yaml', workdir='../chipseq')

    # Create Snakemake subworkflows
    subworkflow rnaseq:
        configfile: rnaseq_config.path
        workdir: rnaseq_config.workdir

    subworkflow chipseq:
        configfile: chipseq_config.path
        workdir: chipseq_config.workdir

    # Here's how to use files created from each workflow:
    rule peaks_at_promoters:
        input:
            de_results=rnaseq(utils.flatten(rnaseq_config.targets['downstream'])),
            peaks=chipseq(utils.flatten(chipseq_config.targets['peaks']))
        output:
            'ma_plot.png'


References workflow and the reference dict
------------------------------------------
the references workflow in
``workflows/references/Snakefile``, when run on its own, builds all references
specified in the config file. This is typically done only when initially
setting up a system that will run workflows on many different references.

Most of the time, this workflow is included (with the ``include:`` directive)
into the other workflows. This way, any reference files that are needed by,
say, the RNA-seq workflow will be created automatically.

The format of the config YAML is designed to be convenient to edit and
maintain, but ends up being awkward to use in practice. So in the ``SeqConfig``
object (see description of this above), we store a more convenient
representation of it that only contains the generated reference files. It is
available as the ``refdict`` attribute.  It converts this, as entered into the
config yaml (and which includes URLs and other info):

.. code-block:: yaml

    references:
      dm6:
        r6-11:
          fasta:
            url: "https://url/to/dm6.fasta"
            indexes:
              - bowtie2
              - hisat2
          gtf:
            url: "https://url/to/gm6.gtf"
            conversions:
              - refflat
        r6-11_transcriptome:
          fasta:
            url: "https://url/to/transcriptome.fa"
            indexes:
              - salmon

to this simplified version where values are filenames:

.. code-block:: python

    {
      'dm6': {
         'r6-11': {
             'fasta': '/data/dm6/r6-11/fasta/dm6_r6-11.fasta',
             'refflat': '/data/dm6/r6-11/gtf/dm6_r6-11.refflat',
             'gtf': '/data/dm6/r6-11/gtf/dm6_r6-11.gtf',
             'chromsizes': '/data/dm6/r6-11/fasta/dm6_r6-11.chromsizes',
             'bowtie2': '/data/dm6/r6-11/bowtie2/dm6_r6-11.1.bt2',
             'hisat2': '/data/dm6/r6-11/hisat2/dm6_r6-11.1.ht2',
             },
         'r6-11_transcriptome': {
             'fasta': '/data/dm6/r6-11_transcriptome/fasta/dm6_r6-11_transcriptome.fasta',
             'chromsizes': '/data/dm6/r6-11_transcriptome/fasta/dm6_r6-11_transcriptome.chromsizes',
             'salmon': '/data/dm6/r6-11_transcriptome/salmon/dm6_r6-11_transcriptome/hash.bin,
             },
        },
    }
