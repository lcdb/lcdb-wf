.. _patterns-and-targets:

Patterns and targets
====================
We use a hybrid approach to specifying input and output patterns for Snakemake
rules. These are used in the RNA-seq and ChIP-seq workflows.

Patterns
--------

**Patterns** are filename-like strings with placeholders in them, like this::

    data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz

Generally you don't need to modify them (though you can). The patterns file is
useful as a guide to the files that are created by the workflow.

- RNA-seq patterns are in ``workflows/rnaseq/config/rnaseq_patterns.yaml``
- ChIP-seq patterns are in ``workflows/chipseq/config/chipseq_patterns.yaml``


Targets
-------
The metadata (sample names, peak-calling runs) configured in the sample table
and config file fills in the patterns to create the targets. It's the
equivalent of a complicated `expand()` call in a standard
Snakefile.

If we had 2 samples, A and B, then filling in the pattern::

    data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz

would result in::

    data/rnaseq_samples/A/A_R1.fastq.gz
    data/rnaseq_samples/B/B_R1.fastq.gz


How patterns and targets are used in Snakefiles
-----------------------------------------------
Briefly:

- Each Snakefile has access to an object, ``c``.
- **Patterns** are accessed via the ``c.patterns`` dictionary. The structure of
  ``c.patterns`` matches that of the patterns file. Patterns still have the
  ``{}`` placeholders as written in that file.
- **Targets** are accessed via the ``c.targets`` dictionary. The structure
  matches that of ``c.patterns``, but the placeholders are filled in based on
  other configuration information such that each single string from patterns
  may turn into a list after being filled in using the `snakemake.expand()`
  function.
- We can collapse arbitrary groups of patterns or targets together into
  a flattened list with ``lib.utils.flatten()``.

Here is an example rule that uses patterns and targets from the ``c`` object:

.. code-block:: python

    rule all:
        input:
            c.targets['cutadapt']

    rule cutadapt:
        input:
            c.patterns['fastq']
        output:
            c.patterns['cutadapt']

        shell:
            'cutadapt {input} -o {output}'

In the example above, ``c.patterns['fastq']`` and ``c.patterns['cutadapt']``
have the following values, configured in ``config/rnaseq_patterns.yaml``:

.. code-block:: python

    c.patterns['fastq']
    # data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz

    c.patterns['cutadapt']
    # data/rnaseq_samples/{sample}/{sample}_R1.cutadapt.fastq.gz


And ``c.targets[['cutadapt']`` might have the following values, after being
filled in with the sample table (see below for details):

.. code-block:: python

    c.targets['cutadapt']
    # data/rnaseq_samples/sample1/sample1_R1.cutadapt.fastq.gz
    # data/rnaseq_samples/sample1/sample2_R1.cutadapt.fastq.gz
    # data/rnaseq_samples/sample1/sample3_R1.cutadapt.fastq.gz
    # data/rnaseq_samples/sample1/sample4_R1.cutadapt.fastq.gz




This has several advantages:

- Patterns can be automatically filled in by the sample table and config file,
  so the workflow is largely controlled by a TSV file and a YAML file.

- Storing filenames outside individual Snakefiles allows us to access them much
  more easily from other Snakefiles or downstream scripts.

- Writing aggregation rules is much easier. For example, instead of lots of
  ``expand()`` calls, we can get all the FastQC output across raw, trimmed, and
  aligned runs with ``flatten(c.targets["fastqc"])``.

- Toggling entire sections of the workflow can be performed by changing
  a single line in the first ``all`` rule.

- We can re-organize the output directories only by editing the patterns file
  -- no need to touch the Snakefile.

- It's easier to understand the output locations of files by looking at the
  patterns file than it is to scroll through a Snakefile.

- Letting the ``c`` objects do the work of filling in patterns allows complex
  work to be abstracted away, resulting in simpler Snakefiles. For example,
  filling in the output BED files across many arbitrary-named configured
  peak-calling runs gets complicated, but since this is handled transparently
  by the ``c`` object, we can do things like
  ``utils.flatten(c.targets['peaks'])`` to get all the BED files for all
  peak-callers and all peak-calling runs.

.. seealso::

    For more details, the code is the authoritative source.

    In particular:

        - :class:`lib.patterns_targets.SeqConfig`
        - :class:`lib.patterns_targets.RNASeqConfig`
        - :class:`lib.patterns_targets.ChIPSeqConfig`

    In addition, the `figures Snakefile
    <https://github.com/lcdb/lcdb-wf/blob/master/workflows/figures/Snakefile>`_
    demonstrates how the ChIP-seq and RNA-seq patterns and targets can be used
    for downstream work.

