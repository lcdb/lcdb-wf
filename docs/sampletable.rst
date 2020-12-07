.. _sampletable:

Sample tables
=============
Sample tables map sample names to files on disk and provide additional
metadata.It is expected to have a header and be tab-delimited. Empty lines and
lines that start with a comment (``#``) are skipped.

For running new experiments, you will need to write your own sample table. For
running experiments uploaded to SRA, **you can use the SRA sample table
as-is**, with the addition of a new column indicating what you would like to
name each sample. This makes it almost trivial to run arbitrary SRA RNA-seq
data sets! For ChIP-seq data from SRA, the additional columns `antibody`,
`label`, and `biological_material` as described below will need to be added,
but often that information is already in the SRA sampletable so the columns
just need to be renamed.

.. _rnaseq-sampletable:

RNA-seq sample table
--------------------
Here is an example minimal sample table for RNA-seq. It only contains sample
IDs for four samples::

    # Minimal RNA-seq sample table
    sample
    c1
    c2
    t1
    t2

In this minimal example, the original FASTQ files are expected to be at the
locations ``data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz``. That pattern
is configured in the ``config/rnaseq_patterns.yaml`` file if you would like to
change it (see :ref:`patterns-and-targets`). Specifically, the workflow will
expect the following files to already exist (paths relative to the Snakefile)::

    # The above sample table expects these files to exist:
    data/rnaseq_samples/c1/c1_R1.fastq.gz
    data/rnaseq_samples/c2/c2_R1.fastq.gz
    data/rnaseq_samples/t1/t1_R1.fastq.gz
    data/rnaseq_samples/t2/t2_R1.fastq.gz

.. _symlinks:

Symlinking FASTQs
~~~~~~~~~~~~~~~~~

To avoid having to copy or symlink files over into the expected directory
structure, we can instead list the original filenames in a column called
``orig_filename`` and they will be automatically symlinked into
``data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz``. That is, the following
sampletable::

    # Example RNA-seq sample table with original filenames are specified
    sample   orig_filename
    c1       /data/c1.fastq.gz
    c2       /data/c2.fastq.gz
    t1       /data/other/t1.fq.gz
    t2       ../raw-data/t2.fq.gz

Will result in the following symlinks::

    data/rnaseq_samples/c1/c1_R1.fastq.gz -> /data/c1.fastq.gz
    data/rnaseq_samples/c2/c2_R1.fastq.gz -> /data/c2.fastq.gz
    data/rnaseq_samples/t1/t1_R1.fastq.gz -> /data/other/t1.fq.gz
    data/rnaseq_samples/t2/t2_R1.fastq.gz -> ../raw-data/t2.fq.gz

Note that `orig_filename` paths in the sampletable are considered *relative to
the Snakefile*.

Additional metadata
~~~~~~~~~~~~~~~~~~~
For RNA-seq, only the first column and optionally the `orig_filename` column
are used directly by the RNA-seq workflow.

However, the sampletable is imported into the ``downstream/rnaseq.Rmd`` file
(see :ref:`downstream` for more info) so it's often useful to include
any metadata in the sampletable so it's all in one place, and you'll get all
that information imported into R.

For example, with this sample table we would be easily able to use a DESeq
model of ``~condition`` since the condition column will be imported into R.

::

    # Example RNA-seq sampletable with "condition" metadata included
    sample   orig_filename          condition
    c1       /data/c1.fastq.gz      control
    c2       /data/c2.fastq.gz      control
    t1       /data/other/t1.fq.gz   treatment
    t2       /data/other/t2.fq.gz   treatment

Paired-end data
~~~~~~~~~~~~~~~
**Paired-end and single-end data may be mixed in the same sampletable.**
A sample is specified as paired-end using a separate column in the sampletable.
That column can either be named `layout` (easiest if you're writing your own
sample table) or `LibraryLayout` (if you're using an SRA sampletable, in which
case you can leave it as-is). An error will be raised if both columns are
provided.

If one of these columns exists, the values of the column are converted to
lowercase. For each sample, if the value is either `pe` or `paired`, the sample
will be considered paired-end. In all other cases the sample will be considered
single-end.

For paired-end samples that will be symlinked, both `orig_filename` and
`orig_filename_R2` must be specified as paths relative to the Snakefile (see
:ref:`symlinks` above). If there is a mix of SE and PE samples, the SE sample
must have an empty entry for `orig_filename_R2` (in the context of the
tab-delimited sampletable, this means two tab characters next to
each other with nothing in between).

.. note::

  If the sample table contains both single- and paired-end samples, the
  `fastq_dump` and `cutadapt` rules will create empty R2 files.

  Once the BAM files are created (after alignment in a single- or paired-end
  fashion as appropriate for the sample), we operate mostly on the BAM.

  After the alignment stage, remaining rules **do not** differentiate between
  single- and paired-end reads. In particular, featureCounts and bamCoverage
  may need different parameters depending on the library layout.

::

    # Example RNA-seq sample table with original filenames are specified,
    # and c1 is a paired-end sample
    sample   orig_filename         orig_filename_R2      layout
    c1       /data/c1_R1.fastq.gz  /data/c1_R2.fastq.gz  PE
    c2       /data/c2.fastq.gz                           SE
    t1       /data/other/t1.fq.gz                        SE
    t2       /data/other/t2.fq.gz                        SE

.. _chipseq-sampletable:

ChIP-seq sample table
---------------------
**Three additional columns are required** for ChIP-seq: ``antibody``,
``biological_material`` and ``label``.


:``antibody``:
    Used for differentiating between input and IP samples. Input samples should
    be listed with an antibody of exactly ``input``.

:``biological_material``:
    Ties together which samples came from the same chromatin. This is how we
    know a particular input sample is the matched control for a particular IP
    sample. This is primarily used in the `fingerprint` rule, where we collect
    all the input BAMs together for performing QC. See the
    `lib.chipseq.merged_input_for_ip` function for the technical details of how
    this is handled.

:``label``:
    Used to tie together technical replicates, and **used to configure the
    ChIP-seq peak-calling runs** (see :ref:`cfg-chipseq`).

    Technical replicates share the same label. If you don't have technical
    replicates, then this column can be a copy of the first column containing
    sample names. Technical replicates will have their BAMs merged together
    and duplicates removed from the merged BAM.

The reason that the ChIP-seq sample table is more complicated than RNA-seq is
because RNA-seq is often analyzed in R, and complicated sample handling (like
summing technical replicates) can be performed very flexibly in R. In contrast,
ChIP-seq peak-callers are command-line tools and frequently only take a single
biological replicate, and so are run as Snakemake rules. As a result, more
complex configuration is required to ensure complex experimental designs are
handled correctly.


Minimal ChIP-seq sample table, no replicates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A minimal ChIP-seq sample table, with no biological replicates, looks like
this::

    # Example minimal ChIP-seq sample table
    sampleid    antibody   biological_material  label          orig_filename
    ip1         gaf        s2cell-1             s2cell-gaf-1   /data/ip1.fastq.gz
    input1      input      s2cell-1             s2cell-input-1 /data/input.fastq.gz

- The input sample is required to have the antibody as "input"
- For an IP, its matched input is the sample with ``antibody == input`` that
  also has the same biological material as the IP. Here, we know `input1` goes
  with `ip1` because they both have the same biological material.


ChIP-seq sample table, biological replicates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is another example, this time with biological replicates::

    # Example ChIP-seq sampletable with biological replicates
    sampleid    antibody   biological_material  label          orig_filename
    ip1         gaf        s2cell-1             s2cell-gaf-1   /data/ip1.fastq.gz
    ip2         gaf        s2cell-2             s2cell-gaf-2   /data/run2/ip3.fastq.gz
    input1      input      s2cell-1             s2cell-input-1 /data/input.fastq.gz
    input2      input      s2cell-2             s2cell-input-2 /data/run2/input2.fastq.gz

- As before, `ip1` and `input1` share the same biological material, indicating
  that `input1` is the matched input for `ip1`.
- The matched input for `ip2` is `input2` because they share the same
  biological material (`s2cell-2`) and `input2` has ``antibody == input``.
- Each sample has a unique label because there are no technical replicates
  here.

ChIP-seq sample table, biological and technical replicates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another example, this time with biological and technical replicates:

::

    # Example ChIP-seq sampletable with bio and tech reps
    sampleid    antibody   biological_material  label          orig_filename
    ip1         gaf        s2cell-1             s2cell-gaf-1   /data/ip1.fastq.gz
    ip1a        gaf        s2cell-1             s2cell-gaf-1   /data/ip2.fastq.gz
    ip2         gaf        s2cell-2             s2cell-gaf-2   /data/run2/ip3.fastq.gz
    input1      input      s2cell-1             s2cell-input-1 /data/input.fastq.gz
    input2      input      s2cell-2             s2cell-input-2 /data/run2/input2.fastq.gz


- `ip1` and `ip1a` are technical replicates because they share the label
  `s2cell-gaf-1`. This is often the case when we need to sequence the same
  sample again for higher depth.

- `ip1` and `ip1a` will be merged into one BAM file named after their common
  label, `s2cell-gaf-1` (described further below). The remaining `ip2`,
  `input1`, and `input2` do not have to be merged with anything, so they will
  be symlinked.

Merging technical replicates for ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In contrast to technical replicates in RNA-seq, where counts can be summed in
R, ChIP-seq is a bit more complicated. The ChIP-seq workflow uses ``samtools
merge`` to merge together the unique, duplicates-removed BAM files from
technical replicates into a single BAM, and then removes the duplicates again
from that merged file.

There is a "merged_techreps" key in ``config/chipseq_patterns.yaml`` which
defines the filenames to which technical replicates will be merged. By default
this pattern is
``data/chipseq_merged/{label}/{label}.cutadapt.unique.nodups.merged.bam``.
After trimming, aligning, removing multimappers, and removing duplicates, tech
reps are merged together. Specifically, these files:

::

    data/chipseq_samples/ip1/ip1.cutadapt.unique.nodups.bam
    data/chipseq_samples/ip1a/ip1a.cutadapt.unique.nodups.bam

get merged and then duplicates removed again from that merged file, resulting
in this file::

    data/chipseq_merged/s2cell-gaf-1/s2cell-gaf-1.cutadapt.unique.nodups.merged.bam

For samples with no technical replicates, only symlinks are performed, so for
example this file::


    data/chipseq_samples/ip2/ip2.cutadapt.unique.nodups.bam

will get symlinked to this file::

    data/chipseq_merged/s2cell-gaf-2/s2cell-gaf-2.cutadapt.unique.nodups.merged.bam

For peak-calling (see :ref:`cfg-chipseq`) and any other downstream analysis,
**the files to use are these merged (or symlinked) BAM files.**
