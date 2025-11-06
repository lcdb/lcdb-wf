
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

Within a workflow directory, the default config file is expected to be at :file:`config/config.yaml`.

Config files, at a minimum, specify which reference FASTA to use (:ref:`reference-config`).

For RNA-seq (:ref:`rnaseq-config`) the config file also specifies a GTF
reference and strandedness of the libraries.

For ChIP-seq (:ref:`chipseq-config`) the config file specifies peak-calling runs.

You can override the default config file location when calling snakemake like
this::

  snakemake --configfile="otherdir/myconfig.yaml" ...

Snakemake will merge the config file(s) given on the command line with the
default config file (:file:`config/config.yaml`).

.. _reference-config:

References section
~~~~~~~~~~~~~~~~~~

This section is just about the references part of the config; see
:ref:`rnaseq-config` and :ref:`chipseq-config` for any additional config for
those workflows.

Using included reference config templates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The repository includes pre-configured reference genome and annotation
templates in :file:`include/reference_config_templates/` for common model
organisms. These templates provide organism name, genome FASTA URL, and
annotation GTF URL (for RNA-seq). They can be used for both ChIP-seq and
RNA-seq to conveniently fill in the references part of the config.

This is the easiest way to configure references. There are two ways to use
these templates:

1. Command-line: Point to the template using ``--configfile`` when calling Snakemake::

    snakemake --configfile=../../include/reference_config_templates/Homo_sapiens/GENCODE.yaml ...

   This merges the template with your default :file:`config/config.yaml`,
   creating new or replacing existing keys.

2. Copy-paste: Copy the contents from a template file into your
   :file:`config/config.yaml` file.

Otherwise, see the next section for customizing the references section.

Configuring references
^^^^^^^^^^^^^^^^^^^^^^

Both RNA-seq and ChIP-seq need a reference fasta configured, like this:

.. code-block:: yaml

  genome:
    url: <URL to gzipped FASTA file>


RNA-seq also needs a GTF annotation configured, which works similarly:

.. code-block:: yaml

  annotation:
    url: <URL to gzipped GTF>


The value of ``url`` can be a file (like
``file:///data/references/Homo_sapiens/gencode.fa.gz``) or any FTP or HTTP URL.

This is useful if you have existing reference files you want to use.

By default, reference files will be downloaded to the :file:`references`
directory within the current workflow. Aligner indexes will be built here as well.

For ChIP-seq, the references directory will look like:

.. code-block:: text

  references/
  ├── genome.fa.gz  # Downloaded FASTA
  ├── bowtie2/  # bowtie2 index
  │   └── genome.*.bt2
  └── genome.chromsizes  # chromsizes from fasta

For RNA-seq, it will look like:

.. code-block:: text

  references/
  ├── bowtie2/  # bowtie2 index for rRNA
  │   └── rrna.*.bt2
  ├── salmon/  # salmon index
  ├── star/   # STAR index
  ├── annotation.gtf.gz  # downloaded GTF
  ├── annotation.refflat # GTF converted to refflat
  ├── annotation.bed12   # GTF converted to bed12
  ├── annotation.mapping.tsv.gz  # TSV of attributes from GTF
  ├── genome.fa.gz  # downloaded FASTA
  ├── genome.fa.fai  # chrom sizes
  ├── rrna.fa.gz  # rRNA sequence for organism from SILVA
  └── transcriptome.fa.gz  # created from genome FASTA and GTF


See :ref:`decisions-references` for a discussion on why it's done this way. You
can control this behavior by using the optional ``references`` entry in the
config file, which will instead look for (and create if needed) the specified
directory. If you do this, keep in mind that each reference directory uses
generic labels like ``genome``, ``annotation``, etc, so using the same
directory for different organisms will cause the files to be overwritten for
the last-run organism. So if you use this approach you should consider putting
your references in directories named after organisms and the versions of
aligners used.



.. _rnaseq-config:

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

In :file:`include/reference_config_templates` you can find configs for common
model organisms. These have both genome and annotation, so you can point
Snakemake to them on the command line. You would still need to specify
strandedness, which can be a config entry in
:file:`config/config.yaml`. Or it could be specified directly on the command
line, like this:

.. code-block:: bash

  snakemake \
    --configfile=../../include/reference_config_templates/Homo_sapiens/GENCODE.yaml \
    --config stranded=fr-firststrand

(in this case a separate :file:`config/config.yaml` would not be needed, as
long as you use the default :file:`config/sampletable.tsv` as your sampletable)

.. _chipseq-config:

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

or don't edit the config and instead modify the command-line call:

.. code-block:: bash

  snakemake --config sampletable="config/mytable.csv" ...

.. _rnaseq-sampletable:

RNA-seq sample table
~~~~~~~~~~~~~~~~~~~~

Here is an example minimal sample table for single-end RNA-seq data. The column
``orig_filename`` is required. Paths are **relative to the Snakefile**.

.. code-block:: text

    # Example RNA-seq sample table
    sample   orig_filename
    c1       /data/c1_R1.fastq.gz
    c2       /data/c2_R1.fastq.gz
    t1       other-data/treatment_1_1.fq.gz
    t2       ../../raw-data/t2_1.fq.gz

For paired-end data, we need to specify the second end of the pair in the
``orig_filename_R2`` column:

.. code-block:: text

    sample   orig_filename                   orig_filename_R2
    c1       /data/c1_R1.fastq.gz            /data/c1_R2.fastq.gz
    c2       /data/c2_R1.fastq.gz            /data/c2_R2.fastq.gz
    t1       other-data/treatment_1_1.fq.gz  other-data/treatment_1_2.fq.gz
    t2       ../../raw-data/t2_1.fq.gz       ../../raw-data/t2_2.fq.gz


The paired-end example above would result in the following symlinks being
created:

.. code-block:: text

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
paired-end sample table we could use:

.. code-block:: text

    sample   group     replicate  orig_filename                   orig_filename_R2
    c1       control   1          /data/c1_R1.fastq.gz            /data/c1_R2.fastq.gz
    c2       control   2          /data/c2_R1.fastq.gz            /data/c2_R2.fastq.gz
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
- You can split up experiments across multiple sampletables if needed
- Descriptive sample names help with interpreting QC and downstream analysis


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
``merged_label`` value will be merged together. Take the following example:

.. code-block:: text

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
like this:

.. code-block:: text

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

Advanced: post-processing reference files
-----------------------------------------

In some cases, reference files may need to be modified after download. This is
becoming increasingly rare thanks to updates from providers like Ensembl and
NCBI, but sometimes files need to be post-processed. For example:

- only a GFF is available, so it needs to be post-processed into GTF format
- extra chromosomes are included that should be removed
- renaming chromosomes (e.g. to match UCSC Genome Browser nomenclature)
- adding transgenic constructs to FASTA and/or GTF
- removing problematic annotations (like trans-splicing events which some tools have issues with)

To handle these situations, a reference file config can take an optional
``postprocess:`` key. This is a string containing a dotted name referring to
a Python function importable by :file:`lib.utils`. For
:file:`lib/postprocess/__init__.py` has many such functions, but you can write
your own.

This is a bit of an advanced topic. See the help and comments in
``lib.utils.download_and_postprocess`` (in the file :file:`lib/utils.py`) for
details; the following attempts to provide enough information and direction
for you to implement your own customizations.

A function used for post-processing must have the signature:

.. code-block:: python

    def name_of_function(tmpfiles, outfile, **kwargs):
        pass

It should expect ``tmpfiles`` to be a list of files that were just downloaded,
and ``outfile`` is the final gzipped file to create.

If the function does not need any kwargs, configure it like this:

.. code-block:: yaml

  genome:
    url: <URL>
    postprocess: "lib.postprocess.name_of_function"

If it needs kwargs, configure it like this:

.. code-block:: yaml

  genome:
    url: <URL>
    postprocess:
      function: "lib.postprocess.name_of_function"
      kwargs:
        kwarg1:
          - list
          - of
          - items
        verbose: true

Note these examples use the genome fasta, but the functionality works for
annotations as well.

If a post-processing function has a keyword argument with starts and ends with
a double underscore (``__``), the config system will assume this is a string
that should be interpreted as a dotted function name and the actual function
will be resolved and passed to the post-processing function.

Here is a complete (and complex) example to illustrate the mechanism. In this
example, we want to include ERCC spike-ins to the reference genome as well as
to the GTF file so we can quantify them. However, only a GFF file is available
for *S. pombe*, so we also need to post-process that into a GTF before
appending ERCC annotations. As another wrinkle, there is no ERCC spike-in GTF,
so we need to create our own from the FASTA file. Here is how this would be
configured:

.. code-block:: yaml

    genome:
      url:
        # S. pombe fasta
        - "ftp://ftp.ensemblgenomes.org/pub/fungi/release-41/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna_sm.toplevel.fa.gz"
        # ERCC fasta
        - "https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA"
      postprocess:
        # See lib/postprocess/ercc.py
        function: "lib.postprocess.ercc.add_fasta_to_genome"

    annotation:
      url:
        # S. pombe GFF, which needs to be converted to GTF
        - "ftp://ftp.ensemblgenomes.org/pub/fungi/release-41/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.41.gff3.gz"

        # ERCC GTF is not available; conversion function needed to convert
        # fasta to GTF
        - "https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA"

      postprocess:
        function: "lib.postprocess.ercc.add_gtf_to_genome"
        kwargs:
          # As per the docs for add_gtf_to_genome(), this function will be
          # applied to all but the last input file. It is specified as a string
          # here, but the config-processing system will resolve this to the
          # actual function and pass that along to add_gtf_to_genome
          __preprocess__: "lib.common.gff2gtf"

The end result is a genomic fasta with ERCC spike-ins added and a GTF version
of Ensembl's GFF file with ERCC spike-ins added as additional annotations.
