``references.snakefile``
========================

When run on its own, this builds all references specified in the config file.
This is typically done when initially setting up a system that will run
workflows on many different references. For smaller-scale use-cases (e.g.,
you'll only be running data from a single organism), you might be better off
not running the workflow directly. This workflow is also included (with the
``include:`` directive) into the other workflows. This way, any reference files
that are needed by, say, the RNA-seq workflow will be created automatically.

See :ref:`references` for more details.

:func:`lib.common.reference_dict`

TODO:

- references dict and how to access
- wrappers
- orig_filename
- R_rnaseq vs top-level env
- rnaseq.Rmd
- test-settings
- trackhub


Features common to all workflows
================================

`patterns` and `targets`
------------------------
At the top of each workflow (snakefile), there is a `patterns` dict that lays
out, all in one place, the output file patterns used throughout the workflow.
Subdictionaries are used for better organization.  Rules can optionally use
values from this dict to use as their input and output patterns. Using the
`lcdblib.snakemake.fill_patterns` function, these patterns are also "rendered"
into a `targets` dictionary (this is basically a recursive `expand()`). So the
`patterns` dictionary has the patterns with wildcards, and `targets` has the
patterns with the wildcards filled in. The `all` rule uses the contents of the
`targets` dict to populate its input, and therefore control what rules make it
into the DAG 

This provides, all in one section of the Snakefile, a fair amount of
customization options. Patterns can be commented out, or targets can be excluded
from the `all` rule to fine-tune which rules will be run.

A nice side-effect of the `patterns` and `targets` design is that aggregation
rules become much easier to write. If all FastQC runs are under a `fastqc` key
in `targets`, they can all be used as input to a rule using
`lcdblib.utils.flatten(targets['fastqc'])`. This is in contrast to, say, writing
input functions.

Here is a small example, first without the patterns and targets:


.. code-block:: python

    rule all:
        input:
            'data/fastqs/{sample}.fastq.fastqc.html',
            'data/multiqc.html'

    rule align:
        input: 'data/fastqs/{sample}.fastq.gz'
        output: 'data/bams/{sample}.bam'
        shell:
            pass

    rule fastqc:
        input: '{prefix}'
        output: '{prefix}.html'
        shell:
            pass

    rule multiqc:
        input:
            expand('data/fastqs/{sample}.fastq.fastqc.html', sample=['1', '2']),
            expand('data/bams/{sample}.bam.fastqc.html', sample=['1', '2']),
        output: 'data/multiqc.html'


With the files hard-coded throughout the snakefile, modifying the locations of
files becomes harder the more complex the file. For example, if we wanted to
change the fastqs location from ``data/fastqs`` to ``data/samples``, we would
have to edit all the rules individually. Furthermore, it's not so easy to get
a high-level picture of the files created by the workflow.

In contrast, here is the same workflow using the patterns and targets. The big
difference is that the filename patterns are set up all in one block at the
beginning so it's easier to see what is created. To change the fastqs location,
we would simply change the pattern for ``fastq:`` and ``fastqc:fastq:``, and
this would propagate to the rules. Last, it is easier to enable/disable
portions of the DAG in the ``all`` rule because the targets are named.

.. code-block:: python

    from lcdblib.snakemake.utils import fill_patterns

    # Dictionary of patterns. Keys and sub-dictionaries are up to the user; the
    # idea is to help with organization.
    patterns = {
        'fastq': 'data/fastqs/{sample}.fastq.gz',
        'mapped': 'data/bams/{sample}.bam',
        'fastqc':
            'fastq': 'data/fastqs/{sample}.fastq.fastqc.html',
            'bam': 'data/bams/{sample}.bam.fastqc.html',
        'multiqc': 'data/multiqc.html',
        'markdups': 'data/bams/{sample}.nodups.bam',
    }

    # Recursively fills out patterns.
    targets = fill_patterns(patterns, sample=['1', '2'])

    # Comment out the "markdups" line to remove that rule from the DAG
    rule all:
        input:
            targets['fastqc'] +
            # targets['markdups'] +
            targets['multiqc']


    # Patterns can be used as input/output
    rule align:
        input: patterns['fastq']
        output: patterns['mapped']
        shell:
            pass

    # In some cases, it makes sense to NOT use the patterns
    rule fastqc:
        input: '{prefix}'
        output: '{prefix}.html'
        shell:
            pass

    # Aggregation rules can use the already-filled-in `targets`
    rule multiqc:
        input: targets['fastqc']
        output: patterns['multiqc']
