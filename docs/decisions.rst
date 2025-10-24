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

featureCounts all-in-one or individually
----------------------------------------

featureCounts can accept a list of BAMs and run everything in one shot, or can
be run once per sample and then manusally aggregated later. Previously, we
provided all BAMs. However, for paired-end BAMs, featureCounts will internally
name sort each BAM before counting. It does this serially. The result is
possibly substantial memory usage and a lot of time. 

One approach could be to temporarily name-sort BAMs in a separate rule,
conditional on paired-end reads, and the featureCounts rule would need to have
conditional input filenames as well. This adds a little bit of complexity for
the benefit of being able to more finely control resource usage. Another
approach would be to run featureCounts independently on each BAM, allosing it
to name-sort independently each one in parallel, and then manually aggregate
the featureCounts output of each.

Since the conditional inclusion of a namesorted rule was straightforward (a
matter of choosing the input file for featureCounts rule), it made the most
sense to run featureCounts once, providing it all samples.

Selection of reference genomes and annotations
----------------------------------------------

Where possible, we select "primary" assemblies -- those with th canonical
chromosomes and unassembled contigs (scaffolds) but NOT haplotypes, alternate
loci, or assembly patches.

`Heng Li's blog post
<https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use>`__ on
the subject is a useful guideline. To summarize, we want to exclude alt contigs
/ haplotypes because they may create multimapping issues, and we want to
include unassembled contigs because excluding them will artificially decrease
alignment percentage.

Since lcdb-wf is intended to be used with arbitrary organisms, the PAR and
mitochondrial sequences mentioned there are not relevant in general.

Ideally, we would have a tool that, given the URLs for raw fastq and gtf,

1. Displays the set of chromosomes
2. Infers if there are any that look like rDNA or mtDNA
3. Ensures the GTF matches the fasta match chromosomes
4. Accepts a template config to assess to process


Annotations
-----------

We use the most comprehensive annotations. For human and mouse, this is the
GENCODE "comprehensive" annotation for the primary assembly, which will include
many more than just protein-coding transcripts. For example, here are the
frequencies of ``transcript_type`` values in GENCODE v19's comprehensive
annotation:

::

  1726632 protein_coding
  214952 nonsense_mediated_decay
  154780 processed_transcript
  135772 retained_intron
   54584 lincRNA
   44207 antisense
   22976 processed_pseudogene
   15313 pseudogene
   11202 unprocessed_pseudogene
    9477 miRNA
    7090 transcribed_unprocessed_pseudogene
    6149 misc_RNA
    5783 snRNA
    4521 snoRNA
    3148 sense_intronic
    1662 polymorphic_pseudogene
    1610 rRNA
    1430 unitary_pseudogene
    1417 sense_overlapping
    1117 IG_V_gene
    1091 transcribed_processed_pseudogene
    1035 non_stop_decay
     755 TR_V_gene
     681 IG_V_pseudogene
     300 TR_J_gene
     185 IG_C_gene
     152 IG_D_gene
     100 3prime_overlapping_ncrna
      99 TR_V_pseudogene
      80 IG_J_gene
      66 Mt_tRNA
      56 TR_C_gene
      36 IG_C_pseudogene
      12 TR_J_pseudogene
      12 TR_D_gene
       9 IG_J_pseudogene
       6 Mt_rRNA
       3 translated_processed_pseudogene

Erring on the side of too many annotations (i.e., using the comprehensive
annotation instead of a curated version) will result in more features, which at
face value might make the FDR adjustment more harsh in DESeq2. But DESeq2's
independent filtering (not even testing those features with so few reads that
they would not reach significance) guards against this.

Zipping/unzipping references
----------------------------

STAR requires uncompressed FASTA and GTF files to build the index. Making
uncompressed files temporary means running the risk of another rule needing
uncompressed to trigger costly STAR alignment. The extra storage cost of
leaving an uncompressed fasta (~3 GB) around is minimal compared to the scale
of all other data, and guards against inadvertently re-running all alignment
jobs.

Test framework
--------------

I had previously thought that the CircleCI tests were annoying to run and
reproduce locally, so the ``tests/lcdb-wf-test`` script was born. Turns out
that got rather complicated, and ended up being just as annoying. In the spirit
of reducing complexity, that test harness script is removed.
