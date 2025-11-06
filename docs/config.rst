
.. _config:


Configuration
=============

Configuration happens in two places:

**Config file:**

- :ref:`rnaseq-config`
- :ref:`chipseq-config`

**Sampletable:**

- :ref:`rnaseq-sampletable`
- :ref:`chipseq-sampletable`


.. _configfiles:

Config file
-----------

Config files, at a minimum, specify which reference FASTA to use (:ref:`reference-config`).

For RNA-seq (:ref:`rnaseq-config`) the config file also specifies strandedness.

For ChIP-seq (:ref:`chipseq-config`) the config file specifies peak-calling runs.

Config files are in YAML format. By default, they are expected to be at
:file:`config/config.yaml`, but you can override from the command line like this::

  snakemake --configfile="otherdir/myconfig.yaml" ...

Snakemake will merge the config file(s) given on the command line with the
default config file (:file:`config/config.yaml`).

.. _reference-config:

Configuring genome fasta (RNA-seq & ChIP-seq)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both RNA-seq and ChIP-seq need a reference fasta configured, like this:

.. code-block:: yaml

  genome:
    url: <URL to gzipped FASTA file>

The value of ``url`` can be a file, like
``file:///data/references/Homo_sapiens/gencode.fa.gz``, or any FTP or HTTP URL.


You could optionally use the included reference configs to fill in the genome
and annotation from the commandline, and Snakemake would be called like this::

  snakemake --configfile=../../include/reference_config_templates/Homo_sapiens/GENCODE.yaml ...

Or you could copy the contents of the reference config templates and paste in
your own :file:`config/config.yaml`.


- url can be file
- postprocessing
- overrides
- included reference configs


RNA-seq config
~~~~~~~~~~~~~~

For RNA-seq, in addition to the genome fasta file described above, you also need:

- ``annotation``, structured similar to ``genome``, which specifies a gzipped
  GTF file. A transcriptome fasta is automatically built from the genome fasta
  and this GTF.
- ``organism`` which will be used to screen ribosomal RNA. Technically, this is
  searching for the string in the SILVA rRNA database's fasta records.
- ``stranded`` of the libraries, which is used for automatically
  configuring strand-specific tools. The options are:
  - ``fr-firststrand`` for dUTP libraries
  - ``fr-secondstrand`` for ligation libraries
  - ``unstranded`` for libraries without strand specificity.

See https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings for more
info on strandedness. If you don't know ahead of time, you can use
``fr-firststrand`` and inspect the results for RSeQC's infer_experiment in the
MultiQC output. Correct the strandedness in the config, and re-run. Only the
jobs affected by strandedness will be re-run.

Here is an example for human:

.. code-block:: yaml

    organism: "Homo sapiens"
    genome:
      url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
    annotation:
      url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"
    stranded: "fr-firststrand"

In :file:`include/reference_configs` you can find configs for common model
organisms. These have both genome and annotation, so you can point Snakemake to
them on the command line. You would still need to specify strandedness, which
can be the only config entry in :file:`config/config.yaml`. Or it could be
specified directly on the command line, like this:

.. code:block:: bash

  snakemake \
    --configfile=../../include/reference_configs/Homo_sapiens/GENCODE.yaml \
    --config stranded=fr-firststrand

(in this case no separate :file:`config/config.yaml` would be needed, as long
as you use the default :file:`config/sampletable.tsv` as your sampletable)


ChIP-seq config
~~~~~~~~~~~~~~~

For ChIP-seq, in addition to the genome fasta file described above, you also
need a peak-calling section if you want to to run peak-calling.

The idea is that the ``peak_calling:`` entry in the config is a list. Each item
in the list is a dictionary with the following keys:

- ``label`` for the peak-calling run. This is intentionally free-form since you
  may want to run the same samples through multiple algorithms or different
  parameters. Output will be in :file:`data/peak_calling/<algorithm>/<label>`.
- ``algorithm``, currently supported options are ``macs``, ``epic2``
- ``ip`` a list of IP samples (or merged tech reps) or their equivalent. For
  ATAC-seq, this is the ATAC samples. For CUT&RUN and Cut&Tag, these are the
  samples with the antibody of interest.
- ``control`` a list of control samples (or merged tech reps). For ATAC-seq,
  leave this empty or exclude entirely. For CUT&RUN and Cut&Tag, this can be
  excluded or IgG samples can be used.
- ``extra`` is a string that is passed verbatim to the peak-caller's
  command-line call. Use this to modify parameters for each peak-calling run.

.. note::

   The values for ``ip`` and ``control`` must match the sampletable's ``label``
   column, see :ref:`chipseq-sampletable` for details. That is, these are the
   names of the merged technical replicates.


Here is an example:

.. code-block:: yaml

    chipseq:
      peak_calling:
        - label: gaf-wingdisc-pooled
          algorithm: macs
          ip:
            - gaf-wingdisc-1
            - gaf-wingdisc-2
          control:
            - input-wingdisc-1
            - input-wingdisc-2
          extra: '--nomodel --extsize 147'

        - label: gaf-wingdisc-pooled-1
          algorithm: epic2
          ip:
            - gaf-wingdisc-1
          control:
            - input-wingdisc-1

        - label: gaf-wingdisc-pooled-no-control
          algorithm: epic2
          ip:
            - gaf-wingdisc-1
            - gaf-wingdisc-2
          control: []
          extra: ''

.. _sampletables:

Sampletable
-----------

Sample tables map sample names to files on disk and provide additional
metadata.

.. note::

    Running data from SRA? You can use the SRA metadata sample table (see
    `example <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP423046>`__)
    as-is, though you may want to add a new column at the beginning with more
    readable sample names. If there are technical replicates, you may need to
    edit the SRA table -- see :ref:`chipseq-sampletable` below.

Sample tables are TSV or CSV (with a respective .tsv or .csv
extension) and have a header. Empty lines and lines that start with a comment
(``#``) are skipped.

**The first column is interpreted as sample name.** FASTQ files specified in the
sampletable will be automatically symlinked to locations named after the sample
name as expected by lcdb-wf (i.e., you don't need to move files anywhere to meet
these expectations).

By default, the Snakefiles look for a file called
:file:`config/sampletable.tsv`. But you can edit the config file or provide
a commandline config option if you want to use something different:

.. code-block:: yaml

    # in the config yaml
    sampletable: "config/mytable.csv"

    # ...remainder of config

or don't edit the config and instead modify the command-line call::

  snakemake --config sampletable="config/mytable.csv" ...

.. _rnaseq-sampletable:

RNA-seq sample table
~~~~~~~~~~~~~~~~~~~~

Here is an example minimal sample table for single-end RNA-seq data. The column
``orig_filename`` is required::

    # Example RNA-seq sample table
    sample   orig_filename
    c1       /data/c1_R1.fastq.gz
    c2       /data/c2_R1.fastq.gz
    t1       other-data/treatment_1_1.fq.gz
    t2       ../../raw-data/t2_1.fq.gz

For paired-end data, we need to specify the second end of the pair in the
``orig_filename_R2`` column::

    sample   orig_filename                   orig_filename_R2
    c1       /data/c1_R1.fastq.gz            /data/c1_R2.fastq.gz
    c2       /data/c2_rR.fastq.gz            /data/c2_R1.fastq.gz
    t1       other-data/treatment_1_1.fq.gz  other-data/treatment_1_2.fq.gz
    t2       ../../raw-data/t2_1.fq.gz       ../../raw-data/t2_2.fq.gz


Relative paths are interpreted relative to the Snakefile
(:file:`workflows/rnaseq`), so the paired-end examplea above would result in the
following symlinks being created::

  data/rnaseq_samples/c1/c1_R1.fastq.gz --> /data/c1_R1.fastq.gz
  data/rnaseq_samples/c1/c1_R2.fastq.gz --> /data/c1_R2.fastq.gz
  data/rnaseq_samples/c2/c2_R1.fastq.gz --> /data/c2_R1.fastq.gz
  data/rnaseq_samples/c2/c2_R2.fastq.gz --> /data/c2_R2.fastq.gz
  data/rnaseq_samples/t1/t1_R1.fastq.gz --> ../../../other-data/treatment_1_1.fq.gz
  data/rnaseq_samples/t1/t1_R2.fastq.gz --> ../../../other-data/treatment_1_2.fq.gz
  data/rnaseq_samples/t2/t2_R1.fastq.gz --> ../../../../../raw-data/t2_1.fq.gz
  data/rnaseq_samples/t2/t2_R2.fastq.gz --> ../../../../../raw-data/t2_2.fq.gz

This sampletable will be read into the downstream differential expression
analysis, so it's a good idea to add lots of metadata here. Here is a final
paired-end sample table we could use::

    sample   group     replicate  orig_filename                   orig_filename_R2
    c1       control   1          /data/c1_R1.fastq.gz            /data/c1_R2.fastq.gz
    c2       control   2          /data/c2_rR.fastq.gz            /data/c2_R1.fastq.gz
    t1       treatment 1          other-data/treatment_1_1.fq.gz  other-data/treatment_1_2.fq.gz
    t2       treatment 2          ../../raw-data/t2_1.fq.gz       ../../raw-data/t2_2.fq.gz

Here are some additional tips:

- If you have technical replicates, list them as separate lines. This lets us inspect
  the QC of the technical replicates independently. The downstream analysis in
  R has a function for summing counts of technical replicates, so the merging
  can be handled later.
- It's helpful to add as many metadata columns as you can, so that the
  sampletable becomes a single source of truth.
- Creating sampletables is by far the most error-prone step -- it's very easy
  to miss changing an R1 to an R2, for example. Double- and triple-check!
- You can split up experiments


.. _chipseq-sampletable:

ChIP-seq sample table
~~~~~~~~~~~~~~~~~~~~~

A ChIP-seq sampletable needs sample names and original filenames that are
symlinked to expected locations. just like the RNA-seq sampletable described
above.

However, if a sample has technical replicates, they have to be specified in the
sampletable for ChIP-seq. This is in contrast to RNA-seq, where we can simply
sum counts of tech reps in R. See :ref:`decisions-techreps` for details.

Use the ``merged_label`` column to control this. Rows with the same
``merged_label`` value will merged together. Take the following example::

  samplename  merged_label  orig_filename
  ip1a        ip1           /data/run1/ip1.fq.gz
  ip1b        ip1           /data/run2/ip1.fq.gz
  ip2                       /data/run1/ip2.fq.gz
  input1                    /data/run1/input1.fq.gz

In this case, we will get individual QC metrics for technical replicates
``ip1a`` and ``ip1b``. Then they will be merged into a single ``ip1`` (the
merged label) BAM file that is ready for peak-calling. The other samples
(``ip2``, ``input1``) do not have technical replicates.

The merging process uses ``samtools merge`` on the tech reps followed by
Picard ``MarkDuplicates``, saving the result in
:file:`data/chipseq_merged/{merged_label}/{merged_label}.cutadapt.unique.nodups.merged.bam`.
If no merging needs to be done (like for ``ip2`` and ``input1`` here), then we
will get a symlink to the respective BAM file. In this way,
:file:`data/chipseq_merged` always has the complete set of BAMs ready for
peak-calling.

The workflow will automatically fill in missing values in ``merged_label`` with
values from the first column. Or, to be explicit, we could write it all out
like this::

  samplename  merged_label  orig_filename
  ip1a        ip1           /data/run1/ip1.fq.gz
  ip1b        ip1           /data/run2/ip1.fq.gz
  ip2         ip2           /data/run1/ip2.fq.gz
  input1      input1        /data/run1/input1.fq.gz


.. note::

  In prior versions of lcdb-wf, the ``merged_label`` column was called just
  ``label``. This is still supported for backward compatibility; if no
  ``merged_label`` column exists then the workflow will make a new
  ``merged_label`` column out of an existing ``label`` column.

**Use** ``merged_label`` **values when configuring peak-calling.** See
:ref:`chipseq-config` for more on the peak-calling section.

With the sampletable above, a peak-calling config section might then look like this:

.. code-block:: yaml

  chipseq:
    peak_calling:
      - label: ip1-macs
        algorithm: macs
        ip: ip1  # Note the use of merged_label here
        input: input1
      - label: ip2-macs
        algorithm: macs
        ip: ip2
        input: input1

.. note::

  In general, you may find it useful to add an antibody column and a chromatin
  prep column to the sampletable so you know which inputs/controls go with
  which IPs.
