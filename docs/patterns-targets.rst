.. _patterns-and-targets:

Patterns and targets
--------------------
We use a hybrid approach to specifying input and output patterns for Snakemake
rules. These are used in the RNA-seq and ChIP-seq workflows.

**Patterns** are defined in

- ``workflows/rnaseq/config/rnaseq_patterns.yaml`` for RNA-seq
- ``workflows/chipseq/config/chipseq_patterns.yaml`` for ChIP-seq

Generally you don't need to modify them (though you can). This file is useful
as a guide to the files that are created by the workflow.

**Targets** are created by filling in those patterns with metadata
(samplenames, peak-calling runs) as configured in the sample table and config
file. It's the equivalent of a complicated `expand()` call in a standard
Snakefile.

Briefly:

- Each Snakefile has access to an object, ``c``
- **Patterns** are accessed via the ``c.patterns`` dictionary. The structure of
  this dictionary matches that of the patterns file. Patterns still have the `{}`
  placeholders as written in that file.
- **Targets** are accessed via the ``c.targets`` dictionary. Their placeholders
  are filled in. The structure matches that of ``c.patterns``, but since the
  placeholders are filled in, the single string from patterns may turn into
  a list after being filled in using the `snakemake.expand()` function.
- We can collapse arbitrary groups of patterns or targets together into
  a flattened list with ``lib.utils.flatten()``.

For example:

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

In the example above, the targets and patterns would have the following values,
as configured by default in the ``rnaseq_patterns.yaml`` file:

.. code-block:: python

    c.targets['cutadapt']
    # data/rnaseq_samples/sample1/sample1_R1.cutadapt.fastq.gz
    # data/rnaseq_samples/sample1/sample2_R1.cutadapt.fastq.gz
    # data/rnaseq_samples/sample1/sample3_R1.cutadapt.fastq.gz
    # data/rnaseq_samples/sample1/sample4_R1.cutadapt.fastq.gz

    c.patterns['fastq']
    # data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz

    c.patterns['cutadapt']
    # data/rnaseq_samples/{sample}/{sample}_R1.cutadapt.fastq.gz


This has several advantages:

- Patterns can be automatically filled in by the sampletable, so the workflow
  is largely controlled by a TSV.

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

- Downstream tasks can access these patterns and targets as well, avoiding the
  need to retype filename patterns and maintain them across multiple workflows.

- Letting the ``c`` objects do the work of filling in patterns allows complex
  work to be abstracted away, resulting in simpler Snakefiles. For example,
  filling in the output BED files across many arbitrary-named configured
  peak-calling runs gets complicated, but since this is handled transparently
  by the ``c`` object, we do things like ``utils.flatten(c.targets['peaks'])``
  to get all the BED files for all peak-callers and all peak-calling runs.

.. seealso::

    For more details, the code is the authoritative source.

    In particular:

        - :class:`lib.patterns_targets.SeqConfig`
        - :class:`lib.patterns_targets.RNASeqConfig`
        - :class:`lib.patterns_targets.ChIPSeqConfig`

    In addition, the `figures Snakefile
    <https://github.com/lcdb/lcdb-wf/blob/master/workflows/figures/Snakefile>`_,
    which demonstrates how the ChIP-seq and RNA-seq patterns and targets can be
    used for downstream work.

