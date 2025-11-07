Decision log
============

This document keeps track of the reasoning behind various architecture decisions.

.. _decisions-references:

References
----------
Here are use-cases we have that are common enough to warrant supporting:

**References should support multiple workflows (ChIP-seq, RNA-seq, etc)**

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
- Re-using indexes is space- and time-efficient in the short term, but experience has
  shown it to be inefficient in time and reproducibility in the long term.
- Keeping everything in the same deployment directory also helps with the
  archiving process.
- We were hesitant to update the references in the central location due to
  being unsure of what was depending on them.
- Overall, here we make the decision that the time and space cost to re-make
  references for each project is worth the gain in simplicity and isolation.

Arguments for and against a separate references workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RNA-seq, ChIP-seq, and the upcoming variant calling all need to do something
with references, including possibly patching them. We have to deal with this
inherent complexity. It initially made sense to put common rules in the
separate references workflow.

However, only a subset of the rules in the references workflow are actually
shared across RNA-seq and ChIP-seq -- currently, only the bowtie2 index
(genome-wide ChIP-seq alignment; rRNA screening for RNA-seq), the fasta rule,
chromsizes, and the generic unzip rule. The other rules in the <v2.0 references
workflow (gtf, mappings, conversion_bed12, conversion_refflat, kallisto_index,
salmon_index, transcriptome_fasta, star_index, rrna) are all unique to RNA-seq.
So the <v2.0 references workflow is actually mostly an RNA-seq-only references
workflow. It would make more sense to have those RNA-seq-specific rules in the
RNA-seq workflow directly.

Furthermore, much of the complexity is handled in the
lib.utils.download_and_postprocess function, rather than in the workflow rules.
This is the function that downloads, figures out what functions to apply for
post-processing, and outputs the prepared file. We already are using the utils
module separately in the ChIP-seq and RNA-seq workflows, so there's no
additional overhead to import it into the Snakefiles. We can use that function
directly.

Last, having a workflow split across two Snakefiles hampers the ability to
understand the complete workflow.

Taken together, it made more sense to eliminate the references workflow
entirely, and port the rules to the respective workflows.


Selection of reference genomes and annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Where possible, we select "primary" assemblies -- those with th canonical
chromosomes and unassembled contigs (scaffolds) but NOT haplotypes, alternate
loci, or assembly patches.

`Heng Li's blog post
<https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use>`__ on
the subject is a useful guideline. To summarize, we want to exclude alt contigs
/ haplotypes because they may create multimapping issues, and we want to
include unassembled contigs because excluding them would artificially decrease
alignment percentage.

Since lcdb-wf is intended to be used with arbitrary organisms, the PAR and
mitochondrial sequences mentioned there are not relevant in general.

Reference genome and annotation sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lcdb-wf has always been organism-agnostic. It would be nice to have a single
source of all genomics data such that we could pass an organism name and get
back the referencs. But even Ensembl and NCBI are not uniform in their support.
Sometimes primary assemblies are available; sometimes primary chromosome fastas
are available but the top-level is actually primary (rat, Ensembl); A GTF might
not be available (pombe, Ensembl); or only a toplevel assembly is available and
we need to remove the haplotypes and alt loci out (hg19, Ensembl).

Reference nomenclature and directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Options considered:

1. ``references`` (top-level of project, shared by all workflows)
2. ``workflows/<workflow>/references`` (workflow-specific)

The possible location ``workflows/references`` is functionally similar to
top-level ``references`` (in a parent directory of individual workflows) but
references is no longer a workflow so it doesn't make sense to have it right in
the ``workflows`` directory. So this was excluded as an option.

Recall that in lcdb-wf <2.0, we have organism and then tag. For example, we
might have configurations available for different human genome assemblies
(hg19, hg38) and in the central location we needed to differentiate between
them (e.g. ``references/human/hg19/``), which we did with tags.

If we assume a single organism per workflow, which seems reasonable and that
the references are workflow-specific, then we don't need any of this.
``workflows/<workflow>/references/genome.fa`` for example should cover it.

This becomes inefficient in the case where there are multiple workflows, all
for the same organism and all the same workflow type. For example, a project
with chipseq and a two different RNA-seq experiments would have three copies of
the genome fasta. However in such cases, manually creating symlinks can get
around this if space is a problem, and I think it's an acceptable workaround
for the benefit of simplified references more generally.

So we might have something like the following:

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

Zipping/unzipping references
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some tools need uncompressed files, others are fine with compressed. For example,
STAR requires uncompressed FASTA and GTF files to build the index, but bowtie2
can use a compressed fasta. gffread nees uncompressed FASTA and GTF to make
a transcriptome fasta.

Previously, anything using a FASTA or GTF would use the uncompressed version,
and the ``unzip`` rule marked the uncompressed output as temporary. The problem
with this was when we wanted to make a change in featureCounts. Since this used
the temp uncompressed GTF file, the ``unzip`` rule needed to run again...but
that would then trigger the STAR rule to rerun, because it too used that temp
file and it was being changed (well, re-created but that's the same to
Snakemake). As a result, we had to spend the time/resource cost to realign
*everything* and all the downstream jobs after alignment, just to run
featureCounts.

Making the featureCounts rule use the compressed GTF avoids this issue.
However, the transcriptome fasta and the STAR index need the uncompressed
references. During testing, there were multiple times when the entire workflow
needed to run because a file marked as temporary was transiently needed. Upon
closer inspection, this was correct behavior, but it happened enough for subtle
reasons that, to avoid future confusion, we keep both compressed and
uncompressed.

Annotations
~~~~~~~~~~~

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
they would not reach significance) guards against this. So we stick with the
comprehensive annotations when available.

.. _decisions-patterns:

Patterns/targets
----------------

Previously, we had a system of "patterns", which were the string filenames with
wildcards, stored in a separate yaml file. The original idea was to provide
flexibility in reorganizing outputs -- if you didn't like where things were
stored (e.g., if you didn't want your files to always have ``{sample}`` in the
basename), then you could edit that one file and everything would be updated.

The config system would fill in the patterns so that you also had a list of the
filled-in targets. This was mildly convenient for aggregation rules like
multiqc that use lots of inputs.

It was also useful in integrative workflows (e.g., making figures using results
from ChIP-seq and RNA-seq workflows), where you could use the single patterns
yaml file as a record of all the files created. This made writing rules and
keeping track of input/output files a bit easier:

.. code-block:: python

  rule downstream_figure:
      input: c.targets["bam"]
      output: "fig1.pdf"

instead of:

.. code-block:: python

  rule downstream_figure:
      input: expand('../rnaseq/data/{sample}/{sample}.cutadapt.markdups.bam', sample=SAMPLES)
      output: "fig1.pdf"

However, over the years, this system proved to be too obscure, hampering
understandibility of the workflows. I'm not aware of anyone making changes to
it to modify the output locations. In fact, over the years we realized that the
consistency of output directory structure across hundreds of projects is
a *major* benefit, so we specifically *don't* want to change output locations.

So in version 2.0, this system has been completely removed, preferring to use
hard-coded filenames and plenty of ``expand()`` calls.

Advanced users can always create their own patterns yaml files to use in
downstream work.

.. _decisions-params:

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

Take the cutadapt rule for example, where we typically would want to include
the adapters but it's not uncommon to add other arguments. Here
we're working with a simplified, single-end version of it:

.. code-block:: python

  rule cutadapt:
      input:
          fastq='{sample}.fastq.gz"
      output:
          fastq='{sample}.cutadapt.fastq.gz'
      threads: 8
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

Notice that we have some shared arguments (``--nextseq-trim``, ``--overlap``,
``--minimum-length``) as well as a PE-specific adapter argument. Converting
this one to params would be something like the following:

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
principle of "everything must go in params" is not useful because it impacts
clarity.

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
why, imagine if we determined that argument inside the ``run:`` block instead
of in params, and then we changed the config file's stranded entry
(``config["stranded"]``). Even though we would want it to re-run (since the
config changed), this rule would NOT re-run because the *code* didn't change --
Snakemake does not *evaluate* the code in a ``run:`` block to determine if it
changed. However, it *does* evaluate the params. So in this case, it's
necessary to keep the strand argument detection in the params to take advantage
of this behavior, and correctly re-run the rule if the config's strand argument
has changed.

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
  call directly, to visually match the equivalent command-line call and to make
  it clear what should be edited.

Lack of sample-specific parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently if we have samples with different library preps that need different
arguments for cutadapt, then they need to be split into two separate workflow
directories and the Snakefiles edited accordingly to have the correct parameters
for rules.

Supporting sample-specific parameters would certainly be possible. But this
would go against the goal of reducing complexity.

For example, we'd need a location to store multiple sets of parameters (probably
in the config file) and a mechanism to retrieve them based on sample names. This
could be an additional column in the sampletable indicating "parameter sets".
Then we could create a lookup table in the config storing the different
parameter sets, with each set containing parameters for all rules. We'd need to
handle default params in case they weren't specified. Then we'd need to have
each rules' ``params:`` directive do the lookup in a sample-specific manner,
which would be a lookup function in :file:`lib/utils.py`.

Again, this would all be possible. But to reduce complexity it is a deliberate
design choice to opt for a simpler approach: use multiple workflow directories
and edit the respective Snakefiles appropriately. In cases where samples need to
be compared or considered together across the workflows, an additional workflow
can be introduced to aggregate their output.

featureCounts all-in-one or individually
----------------------------------------

featureCounts can accept a list of BAMs and run everything in one shot, or can
be run once per sample and then those outputs can be aggregated later. Previously, we
provided all BAMs to a single all-in-one call of featureCounts. However, for
paired-end BAMs, featureCounts will internally name sort each BAM before
counting. It does this serially. The result is possibly substantial memory
usage and a lot of time.

One approach could be to temporarily name-sort BAMs in a separate rule,
conditional on paired-end reads, and the featureCounts rule would need to have
conditional input filenames as well. This adds a little bit of complexity for
the benefit of being able to more finely control resource usage. Another
approach would be to run featureCounts independently on each BAM, allowing it
to name-sort independently each one (which would happen in parallel jobs
managed by Snakemake), and then manually aggregate the featureCounts output of
each.

Turns out the conditional inclusion of a namesorted rule was straightforward (a
matter of choosing the input file for featureCounts rule), it made the most
sense to run featureCounts once, providing it all samples, and having it use
the temporarily name-sorted BAMs as input for paired-end experiments.


.. _decisions-testframework:

Test framework
--------------

I had previously thought that the CircleCI tests were annoying to run and
reproduce locally, so the ``tests/lcdb-wf-test`` script was born. Turns out
that got rather complicated, and ended up being just as annoying. In the spirit
of reducing complexity, that test harness script is removed. In part, the new
reference config simplification allows control over configs from the
commandline, reducing the need to handle that from the test script.

rRNA
----
Assessing ribosomal RNA contamination is an important QC step. Different
annotation sources have different ways of indicating ribosomal RNA. For example,
Ensembl GTF files typically have "trancript_biotype" attributes on transcript
featuretypes and "gene_biotype" attributes on gene features, depending on
version (older versions have separate rRNA featuretypes). FlyBase uses separate
rRNA feature types. Dictyostelium does not have anything in the GTF. PomBase
uses the "biotype" attribute.

One way of handling this is to have post-processing steps that extract the rRNA
features from a GTF (probably defaulting to assuming an Ensembl-like
"gene_biotype" attribute) and convert them to `IntervalList format
<https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html>`__
to pass to Picard CollectRnaSeqMetrics.

Another way is to bypass the GTF altogether and align to rRNA directly, which is
what we have historically done here. Previously, the reference configs would all
need an rRNA entry that basically did the same thing for each organism, since
every model organism we've worked with is in the SILVA database. It would
download the full SILVA fasta (for large and small subunits), grep out the
records for our species of interest, and build a bowtie2 index out of that. That
means this method is more general, and arguably more complete, but has its own
complexity: we need to download and filter the fasta, build the bowtie2 index,
and aggregate the results into a MultiQC module.

In the 2.0 refactor, rRNA fasta creation now only needs an organism name and the
Snakefile does what was always in the references config, which is to use the
post-process mechanism to filter the fasta.



Aligners
--------

Previously, HISAT2 and STAR were both supported; salmon and kallisto were both
supported. This created additional complexity in the references workflow and in
the configs. Now, we're just using STAR and salmon (for RNA-seq) and bowtie2 for
ChIP-seq.

Aligners don't seem to make that much of a difference, and officially
supporting just one (plus a psueodaligner for RNA-seq) makes the workflows and
config simpler.


.. _decisions-sample-specific-params:


PEP support
-----------

Support for `Portable Encapsulated Projects
<http://pep.databio.org/spec/specification/>`__ is built into Snakemake. Using
a combination of PEP config files, sample tables, and subsample tables, it is
possible to set up the workflows to use PEP in such a way that it can be
backwards-compatible with prior lcdb-wf versions. Specifically, by providing
TSV sampletables, forcing a sample column name, and populating the table with
subsamples. It would be convenient to offload the complexity of handling
technical replicate configuration to a third-party package.

However, getting technical replicates to work correctly proved to be tricky,
due to the way they come in as lists in the resulting dataframe with PEP. While
it would be possible to fix this, some initial experimentation with this
suggested that it would actually be more complex to do that, so deferring to
another package did not result in a net gain in convenience or in complexity
reduction.

PEP configs are not ruled out completely, but we might need a rewiring and
possible rewriting of the ChIP-seq (and possibly RNA-seq) workflows to fully
support PEP subsamples. I don't consider that effort to be worth it right now,
especially because the current config system already supports technical
replicates.

.. _decisions-techreps:

Technical replicates
--------------------
In practice, it's not uncommon for something to go wrong in library prep or
sequencing such that it makes sense to re-do a library. Typically, if it's just
resequencing the same library (perhaps after rebalancing the multiplexing), we
consider that a technical replicate.

The conventional method for handling technical replicates in RNA-seq is to sum
the counts. That is, we take the Salmon or featureCounts files, where technical
replicates are quantified separately, and sum them after import into R. This
allows us to check QC on individual tech reps e.g. to see if they worked. If we
merged at an early stage (like cocatting the FASTQs), then we would not be able
to check QC separately.

For ChIP-seq, the conventional method is to merge BAM files. However, we still
want to keep observability of individual technical replicates where possible,
which includes inspecting duplicates. However, when we merge BAMs of technical
replicates that each had duplicates removed, it's possible that we're
introducing additional duplicates. So we do another round of duplicate removal
after merging.

The end result of all of this is that we get MultiQC output for all of the
technical replicates separately. For ChIP-seq, the post-merging files are
bigWigs and merged-and-deduped BAMs. Currently these do not have separate
entries in MultiQC.

Removing built-in support for plotFingerprint
---------------------------------------------

deepTools' `plotFingerprint
<https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html>`__
needs matched input to each antibody. Previously, we configured this in the
sampletable with the combination of "biological_material" and "antibody"
columns. Samples with exactly "input" as the antibody were the matched control
for the non-input samples with the same biological material.

This ended up being a little complicated because "biological material" is easily
confused with "biological replicate". And now with common CUT&RUN and Cut&Tag
assays that use IgG as control, "IgG" and "control" should probably be aliases
for "input".

It turns out the "biological_material" column was only ever used for the
plotFingerprint rule. It introduced complexity (in code, configuration,
documentation, and user support) for a single rule. In addition, in practice we
ended up visualizing the bigWigs rather than relying exclusively on the
plotFingerprint metrics. So to reduce complexity, plotFingerpring support is
being removed.

Clearer ChIP-seq config
-----------------------

For "label", it was not clear that it was the merged name. And even if there
were no technical replicates in an experiment, it still needed to be filled out
with copies of the sample name.

Now, ``merged_label`` is an alias for ``label``. If the column is missing
entirely, or if the value is empty for a row, then the samplename will be used
automatically.

Removal of autobump
-------------------
For several versions, resources were wrapped with the ``autobump()`` function,
which would automatically retry jobs with more resources if they failed. Turns
out this wasn't as helpful as expected, because errors (like syntax errors or
other mistakes) ended up being a lot more frequent than exceeding resources.
This resulted in escalating resource allocations and longer run time with no
need. So the autobump was removed.

Cleanup of lib/utils.py
-----------------------
We had accumulated a lot of useful functions over time, but things have changed
enough that they haven't been used. To avoid clutter and additional maintenance
burden in supporting otherwise unused code, these functions were removed.


