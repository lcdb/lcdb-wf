.. _sampletable:

Sampletable
===========

The sample table maps sample names to files on disk and provides additional
metadata. It is expected to have a header and be tab-delimited. Empty lines and
lines that start with a comment (``#``) are skipped.

The first column (column name is not checked) should be the sample name. These
sample names will be used to fill in the patterns to create targets (see
:ref:`patterns-and-targets`).

If there is a column in the sample table called `orig_filename`, that column
should contain the path to the raw gzipped FASTQ file for that sample. The path
is relative to the Snakefile. A symlink will be created from this file to
``data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz`` for RNA-seq and
``data/chipseq_samples/{sample}/{sample}_R1.fastq.gz`` for ChIP-seq. This is
configured by the the `fastq` field in the patterns file (see
:ref:`patterns-and-targets` for more on this).

If the `orig_filename` column does not exist, then the gzipped FASTQ files
are expected to already be at
``data/rnaseq_samples/{sample}/{sample}_R1.fastq.gz`` for RNA-seq or
``data/chipseq_samples/{sample}/{sample}_R1.fastq.gz`` for ChIP-seq.

RNA-seq-specific columns
------------------------
For RNA-seq, only the first column and optionally the `orig_filename` column
are used directly by the RNA-seq workflow.

However, the sampletable is imported into the ``downstream/rnaseq.Rmd`` file
(see :ref:`downstream-rnaseq` for more info) so it's often useful to include
any metadata in the sampletable so it's all in one place, and you'll get all
that information imported into R.

For example, with this sample table we would be easily able to use a DESeq
model of ``~condition`` since the condition column will be imported into R.

::

    # Example RNA-seq sampletable
    sample   orig_filename          condition
    c1       /data/c1.fastq.gz      control
    c2       /data/c2.fastq.gz      control
    t1       /data/other/t1.fq.gz   treatment
    t2       /data/other/t2.fq.gz   treatment

.. _chipseq-specific-columns:

ChIP-seq-specific columns
-------------------------

**Three additional columns are required** for ChIP-seq: ``antibody``,
``biological_material`` and ``label``.

``antibody`` is mostly for differentiating between input and IP samples. Input
samples should be listed with an antibody of exactly ``input``.

``biological_material`` is used to tie together which samples came from the
same chromatin. This is how we know a particular input sample is the matched
control for a particular IP sample. This is primarily used in the `fingerprint`
rule, where we collect all the input BAMs together for performing QC. See the
`lib.chipseq.merged_input_for_ip` function for the technical details of how
this is handled.

``label`` is used to tie together technical replicates. Technical replicates
share the same label If you don't have technical replicates, then this column
can be a copy of the first column containing sample names.  Technical
replicates will have their BAMs merged together and duplicates removed from the
merged BAM.

Here is an example where the first biological replicate, `ip1`, had to be
resequenced. It may help to recognize that unique values in the `sampleid`
column correspond to individual FASTQ files.

::

    # Example ChIP-seq sampletable
    sampleid    antibody   biological_material  label          orig_filename
    ip1         gaf        s2cell-1             s2cell-gaf-1   /data/ip1.fastq.gz
    ip2         gaf        s2cell-1             s2cell-gaf-1   /data/ip2.fastq.gz
    ip3         gaf        s2cell-2             s2cell-gaf-2   /data/run2/ip3.fastq.gz
    input1      input      s2cell-1             s2cell-input-1 /data/input.fastq.gz
    input3      input      s2cell-2             s2cell-input-2 /data/run2/input2.fastq.gz

-  `ip1` and `ip2` are technical replicates because they share the label
  `s2cell-gaf-1`
- Their matched input sample is `input1`, because all of them came from the
  same biological material, `s2cell-1`. There is no technical replicate for
  `input1`, so it is the only sample with a label `s2cell-input-1`.
- `ip3` is the same antibody as `ip1` and `ip2`, but came from different
  biological material. It is a biological replicate. It is the only technical
  replicate for this biological replicate, so it is the only sample with
  a label of `s2cell-gaf-2`.
- The matched input sample for `ip3` is `input3` because it came from the same
  biological material. It is the only technical replicate so it is the only
  sample with a label of `s2cell-input-2`.
