Decision log
============

This document keeps track of the reasoning behind various architecture decisions.

References
----------
Here are use-cases we have that are common enough to warrant supporting:

**References should support multiple workflows (ChIP-seq, RNA-seq, etc)**

- This implies that the means the references dir should be in the ``workflows``
  directory or above.
- For example, this may mean a STAR index for RNA-seq, a bowtie2 index for rRNA
  contamination, and another bowtie2 index for ChIP-seq.

**References should support different organisms in different workflows. There
should be only one organism per workflow though.**

- For example, ``workflows/mouse-rnaseq`` and ``workflows/human-rnaseq`` should
  be supported in the same project.


**References should be re-created for each project.**

- Historically we had a central location for the references (shared by multiple
  deployments of lcdb-wf over the years) but we got conflicts where one
  deployment's aligner version was more recent, causing errors when using the
  index for an older version.
- To keep using this, we'd need to version indexes based on aligner version.
- However, when writing up methods for a paper we need to be able to trace
  back what commands were run to generate the reference, including additional
  patching that may have taken place (as is supported by the references
  workflow).
- Re-using indexes is space- and time-efficient in the short term, but has
  shown to be inefficient in time and reproducibility in the long term.
- Keeping everything in the same deployment directory also helps with the
  archiving process. 
- We were hesitant to update the references in the central location due to
  being unsure of what was depending on them.
- Overall, making the decision that the time and space cost to re-make
  references for each project is worth the gain in simplicity and isolation.

Reference nomenclature and directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Options considered:

1. ``references`` (top-level of project, shared by all workflows)
2. ``workflows/<workflow>/references`` (workflow-specific)

The location ``workflows/references`` is functionally similar to top-level
``references`` (in a parent directory of individual workflows) but references
is no longer a workflow so it doesn't make sense to have it right in the
``workflows`` directory.

Recall that in lcdb-wf <2.0, we have organism and then tag. For example, we
might have configurations available for different human genome assemblies
(hg19, hg38) and in the central location we needed to differentiate between
them (e.g. ``references/human/hg19/``).

If we assume a single organism per workflow, and that the references are
workflow-specific, then we don't need any of this.
``workflows/<workflow>/references/genome.fa`` for example should cover it.

This becomes inefficient in the case where there are multiple workflows, all
for the same organism and all the same workflow type. However in such cases,
manually creating symlinks can get around this, and I think it's an acceptable
workaround for the benefit of simplified references more generally.

::

  workflows/rnaseq/references
    genome.fasta
    genome.chromsizes
    rrna.fasta
    annotation.gtf
    annotation.bed12
    annotation.refflat
    transcriptome.fasta
    star/
      genome.fasta <symlink to ../genome.fasta>
      <star files>
    bowtie2/
      rrna.fasta <symlink to ../rrna.fasta>
      <bowtie2 files>
    salmon/
      transcriptome.fasta <symlink to ../transcriptome.fasta>
      <salmon files>

For ChIP-seq:

::

  workflows/chipseq/references
    genome.fasta
    genome.chromsizes
    bowtie2/
      genome.fasta <symlink to ../genome.fasta>
      <bowtie2 files>


Params
------
The ``params:`` directive allows `non-file parameters for rules
<https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#non-file-parameters-for-rules>`__.
Much (perhaps all?) of what can be done in a ``params:`` directive can also be
done in the body of ``run:`` block. On one hand, it can be nice to have a plain
string ``shell:`` block, and put the complexity in the params. But on the other
hand, sometimes it is harder to follow what's happening in params than it would
be in Python in a ``run:`` block.

This section talks about when and why we use params in lcdb-wf.

One of the nice things sbout Snakemake is that the rules (in ``shell:`` blocks)
can be quite close to the equivalent command-line call. Since rules in these
Snakefiles are intended to be edited, it makes sense to keep them as close to
the command-line as is reasonable.

Take the cutadapt rule, for example, where we typically would want to include
the adapters in the call, but it's not uncommon to add other arguments. Here
we're working with a simplified, single-end version of it:

.. code-block:: python

  rule cutadapt:
      input:
          fastq='{sample}.fastq.gz"
      output:
          fastq='{sample}.cutadapt.fastq.gz'
      threads:
          8
      shell:
          "cutadapt "
          "-o {output[0]} "
          "-j {threads} "
          "--nextseq-trim 20 "
          "--overlap 6 "
          "--minimum-length 25 "
          "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
          "{input.fastq[0]} "
          "&> {log}"


Here's an extreme way of adding params where we pull out each argument into
a separate params item. This isn't very flexible and has lots of repetition, so
we probably don't want this::

.. code-block:: python

  rule cutadapt:
      input:
          '{sample}.fastq.gz"
      output:
          '{sample}.cutadapt.fastq.gz'
      threads:
          8
      params:
          nextseq_trim="--nextseq-trim 20",
          overlap="--overlap 6",
          minimum_length=25,
          a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
      shell:
          "cutadapt "
          "-o {output} "
          "-j {threads} "
          "{params.nextseq_trim} "
          "{params.overlap} "
          "{params.minimum_length} "
          "{params.a} "
          "{input} "
          "&> {log}"

But we could add the arguments to be a single "extra" string and store that
in params, like this:

.. code-block:: python

  rule cutadapt:
      input:
          '{sample}.fastq.gz"
      output:
          '{sample}.cutadapt.fastq.gz'
      threads:
          8
      params:
          extra=(
              "--nextseq-trim 20 "
              "--overlap 6 "
              "--minimum-length 25 "
              "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
          )
      shell:
          "cutadapt "
          "-o {output} "
          "-j {threads} "
          "{params.extra} "
          "{input} "
          "&> {log}"

One thing that's nice about this is that the "changeable things" are visually in
a different location. When running Snakemake with `-p` then the params will be
filled in to make one long string, which we could use for debugging.

But we want to support single- and paired-end reads, and the arguments to
cutadapt depend on that. Here's the actual rule:

.. code-block:: python

  rule cutadapt: input:
          fastq=expand("data/rnaseq_samples/{{sample}}/{{sample}}_R{n}.fastq.gz", n=n),
      output:
          fastq=expand(
              "data/rnaseq_samples/{{sample}}/{{sample}}_R{n}.cutadapt.fastq.gz", n=n
          ),
      log:
          "data/rnaseq_samples/{sample}/{sample}_cutadapt.fastq.gz.log",
      threads: 6
      resources:
          mem="2g",
          runtime="2h",
      run:
          if is_paired:
              shell(
                  "cutadapt "
                  "-o {output[0]} "
                  "-p {output[1]} "
                  "-j {threads} "
                  "--nextseq-trim 20 "
                  "--overlap 6 "
                  "--minimum-length 25 "
                  "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
                  "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
                  "{input.fastq[0]} "
                  "{input.fastq[1]} "
                  "&> {log}"
              )
          else:
              shell(
                  "cutadapt "
                  "-o {output[0]} "
                  "-j {threads} "
                  "--nextseq-trim 20 "
                  "--overlap 6 "
                  "--minimum-length 25 "
                  "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
                  "{input.fastq[0]} "
                  "&> {log}"
              )

Notice that we have some shared arguments as well as a PE-specific adapter
argument. Converting this one to params would be something like the following:

.. code-block:: python

  rule cutadapt: input:
          fastq=expand("data/rnaseq_samples/{{sample}}/{{sample}}_R{n}.fastq.gz", n=n),
      output:
          fastq=expand(
              "data/rnaseq_samples/{{sample}}/{{sample}}_R{n}.cutadapt.fastq.gz", n=n
          ),
      log:
          "data/rnaseq_samples/{sample}/{sample}_cutadapt.fastq.gz.log",
      threads: 6
      resources:
          mem="2g",
          runtime="2h",
      params:
          shared=(
               "--nextseq-trim 20 "
               "--overlap 6 "
               "--minimum-length 25 "
               "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
          ),
          se_pe_specific=(
                 "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
          ) if is_paired else ""
      run:
          if is_paired:
              shell(
                  "cutadapt "
                  "-o {output[0]} "
                  "-p {output[1]} "
                  "-j {threads} "
                  "{params.shared} "
                  "{params.se_pe_specific} "
                  "{input.fastq[0]} "
                  "{input.fastq[1]} "
                  "&> {log}"
              )
          else:
              shell(
                  "cutadapt "
                  "-o {output[0]} "
                  "-j {threads} "
                  "{params.shared} "
                  "{params.se_pe_specific} "
                  "{input.fastq[0]} "
                  "&> {log}"
              )

Note in this case we need to provide ``-o`` and ``-p`` arguments
separately for paired-end. So we still need to have the ``if is_paired`` clause
in the body of the rule. This one could be a little bit confusing with the
``se_pe_specific`` clause, but otherwise it supports both SE and PE.

What if we split that out into params as well, so that everything SE or PE
specific is handled there?

.. code-block:: python

  rule cutadapt:
      input:
          fastq=expand(
              "data/chipseq_samples/{sample}/{sample}_R{n}.fastq.gz",
              n=n, allow_missing=True),
      output:
          fastq=expand(
              "data/chipseq_samples/{sample}/{sample}_R{n}.cutadapt.fastq.gz",
              n=n, allow_missing=True),
      log:
          "data/chipseq_samples/{sample}/{sample}_cutadapt.fastq.gz.log",
      threads: 6
      resources:
          mem="2g",
          runtime="2h",
      params:
          extra=(
              "--nextseq-trim 20 "
              "--overlap 6 "
              "--minimum-length 25 "
              "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
          ),
          se_pe_specific=(
              "-o {output[0]} "
              "-p {output[1]} "
              "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
              "{input.fastq[0]} "
              "{input.fastq[1]} "
              if not is_paired else
              "{input.fastq[0]} "
              "-o {output[0]} "
          )
      shell:
          "cutadapt "
          "-j {threads} "
          "{params.se_pe_specific} "
          "{params.extra} "
          "&> {log}"

Now it becomes a little harder to understand what's going on, and we may have
gone too far in pulling everything out into params. So maybe an absolute
principle of "everything in params" is not useful.

Let's take another example, the featureCounts rule for RNA-seq:

.. code-block:: python

  rule featurecounts:
      input:
          annotation=rules.gtf.output,
          bam=rules.markduplicates.output.bam,
      output:
          "data/rnaseq_samples/{sample}/{sample}_featurecounts.txt",
      log:
          "data/rnaseq_samples/{sample}/{sample}_featurecounts.txt.log",
      threads: 8
      resources:
          mem="16g",
          runtime="2h",
      params:
          strand_arg={
              "unstranded": "-s0",
              "fr-firststrand": "-s2",
              "fr-secondstrand": "-s1",
          }[config["stranded"]],
          se_pe_specific=(
            "-p --countReadPairs" if is_paired
            else ""
          ),
          extra="",
      run:
          shell(
              "featureCounts "
              "{params.strand_arg} "
              "{params.se_pe_specific} "
              "{params.extra} "
              "-T {threads} "
              "-a {input.annotation} "
              "-o {output} "
              "{input.bam} "
              "&> {log}"
          )

Here, it is important to have ``strand_arg`` be in the params. To understand
why, imagine if instead we determined that argument inside the ``run:`` block,
and then we changed the config file's stranded entry (``config["stranded"]``).
Then this rule would NOT re-run because the code didn't change -- Snakemake
does not *evaluate* the code in a ``run:`` block to determine if it changed.
However, it *does* evaluate the params. So in this case, it's necessary to keep
the strand argument detection in the params to take advantage of this behavior,
and correctly re-run the rule if the config's strand argument has changed.

Next, we would want to decide whether *all* arguments should go in ``params:``.
In this case, since we're sort of forced to split out ``strand_arg``, we might
as well split everything out.

In the end we have these observations:

- strand-specific arguments *must* be in ``params:``
- some tools have SE/PE-specific arguments. These need an ``if`` clause
  *somewhere*, whether in a ``run:`` block or in ``params:``
- understandability and configuration flexibility are important goals of lcdb-wf
- factoring out *everything* into params weakens understandibility


Guidelines:

- Stranded arguments must be in params
- SE/PE arguments should be handled inside a ``run:`` block
- Any other arguments should be written in a  ``shell:`` block or a ``shell()``
  call directly, to visually match the equivalent command-line call

Arguments for and against a separate references workflow
--------------------------------------------------------

RNA-seq, ChIP-seq, and the upcoming variant calling all need to do something
with references, including possibly patching them. So we have to deal with this
inherent complexity. It initially made sense to put such common rules in the
separate references workflow.

However, only a subset of the rules in the references workflow are actually
shared across RNA-seq and ChIP-seq -- currently, only the bowtie2 index
(genome-wide ChIP-seq alignment; rRNA screening for RNA-seq), the fasta rule,
chromsizes, and the generic unzip rule. The others (gtf, mappings,
conversion_bed12, conversion_refflat, kallisto_index, salmon_index,
transcriptome_fasta, star_index, rrna) are all unique to RNA-seq. So the
current references workflow is actually mostly an RNA-seq-only references
workflow.

Furthermore, much of the complexity is handled in the
lib.utils.download_and_postprocess function, rather than in the workflow rules.
We already are using the utils module separately in the ChIP-seq and RNA-seq
workflows, so there's no additional overhead to import it.

Last, having a workflow split across two Snakefiles hampers the ability to
understand the complete workflow.

Taken together, it made more sense to eliminate the references workflow
entirely, and port the rules to the respective workflows.
