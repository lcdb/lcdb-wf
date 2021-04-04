
.. _references-config:

References config
=================

The references section defines which genomes, transcriptomes, and annotations
to use. It supports arbitrarily many species and assemblies, and supports
customizing references for a particular project. For example, there are tools
and examples for adding ERCC spike-in controls to references, or adjusting
chromosome nomenclature, or removing problematic entries from GTF files, and so
on.

Another advantage is that this makes it easy to add multiple genomes to the
screening step (which uses fastq_screen).

Specifying a references directory
---------------------------------

See :ref:`cfg-references-dir` for more information on where the references are
built (as well as how and why to adjust this).

Including existing reference configs
------------------------------------
We provide a number of pre-configured reference configs for common model
organisms. If you just want to use some references configs that work, then put
this in your config file:

.. code-block:: yaml

    include_references:
      - '../../include/references_configs'

This will populate the config with the contents of all the files contained in
the ``include/references_configs`` directory. Any paths provided under the
``include_references`` key are relative to the Snakefile using the config.
**Note that you will still need to inspect the contents** of those files to
decide which organsim and tag you want to use for your particular experiment
(see :ref:`cfg-organism` and :ref:`cfg-aligner` for more on these fields). For
example, if you are working with human RNA-seq data, and you use the above
``include_references``, you may want this in your config:

.. code-block:: yaml

    organism: 'human'
    aligner:
      tag: 'gencode-v25'
      index: 'hisat2'
    salmon: 'gencode-v25'

The reason for using ``gencode-v25`` is because that tag is configured for the
``human`` key in ``../../include/references_configs/Homo_sapiens.yaml``.

You can provide entire directories of reference configs, a single file, or use
the references section below. The prioritization works like this:

- an organism can show up in multiple configs; if a tag exists for an organism
  in more than one config, higher-priority configs will overwrite the contents
  of the tag.
- directories have lowest priority; when multiple directories are specified the
  last one has priority
- files have priority over directories; when multiple files are specified the
  last one has priority
- the ``references:``` section always has priority over anything in
  ``include_references:``.

The remainder of this section of the documentation explains how to customize
the references, to add your own or modify the existing examples.

Overview
--------
The references workflow is based on the idea that while each genome's source
files may differ, they can usually be modified to a uniform format. For
example, reference files (FASTA, GTF) may come from different providers
(Ensembl, FlyBase, UCSC, etc) and have slightly different formatting (strange
headers, one big file or a tarball of individual chromosomes, etc), once they
are well-formatted they can be used to create a hisat2 index, a bowtie2 index,
a list of genes, intergenic regions, and so on without any further
customization.

The challenging part is the "well-formatted" part. To solve this, the config
file and references building system allows a very flexible specification of how
to modify references via a plugin architecture. It works something like this:

- Each key in the references section refers to an **organism**.
- An organism has one or more **tags**.
- Each **tag** has a FASTA file and/or a GTF file associated with it.
- Each FASTA or GTF specifies one or more URIs from which to download the raw
  file(s). These can be `ftp://`, `http://`, `https://`, or `file://` URIs.
- An optional **postprocess** key specifies the import path to a Python module.
  This is the primary hook for customization, and is described in more detail
  below.
- For FASTA files one or more **indexes** are requested
- For GTF files, zero or more **conversions** are requested.

.. note::

    If using a ``file://`` URI, it needs to be gzipped.

It's probably easiest to show an example config and then describe what's
happening.

Example references config
-------------------------

The following example configures the workflow to:

- download a fasta file from the GENCODE project for the human genome and build
  a hisat2 and bowtie2 index
- download the corresponding GTF file from GENCODE, strip off the dotted
  version numbers from Ensembl gene and transcript IDs, and create a refFlat
  format file from it
- download the SILVA rRNA database and keep only the ribosomal RNA sequence
  corresponding to *Homo sapiens*

This example contains sufficient real-world complexity to illustrate the
flexibility afforded by the references workflow. It is heavily commented for
illustration.

.. code-block:: yaml

    # EXAMPLE REFERENCES CONFIG SECTION

    # This configures the directory in which the prepared references will be
    # saved (see below for directory structure). If you already have reference
    # files saved in the lcdb-wf structure, point this to that directory to
    # avoid rebuilding a fresh set of references:

    references_dir: 'data/references'

    # One of the organisms configured below. We are only configuring a single one
    # so "human" is our only option here:

    organism: 'human'

    # Here we specify which tag under "human" to use for aligning, as well as
    # which index we'll be using. This example is RNA-seq, so we'll use HISAT2:

    aligner:
      tag: 'gencode-v25'
      index: 'hisat2'

    # Top-level section for references:

    references:

      # Label for this organism or species:

      human:

        # "gencode-v25" is our tag to describe this particular FASTA and GTF
        # we're preparing:

        gencode-v25:

          # This block will define how to get and postprocess a FASTA file:

          fasta:

            # URL to download:

            url: 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_organism.genome.fa.gz'

            # We can optionally build indexes for various aligners:

            indexes:
              - 'hisat2'
              - 'bowtie2'

          # This next block will define how to get and postprocess a GTF file.
          # The coordinates of the GTF file correspond to the
          # coordinates in the fasta defined above, so we're putting it under
          # the same tag. This is not required; we could also put it under
          # separate tag (perhaps called "gencode-v25-annotations")

          gtf:
            url: 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz'

            # The GENCODE annotations include the dotted Ensembl versions in
            # the gene IDs. The following function, strip_ensembl_version, is
            # defined in lib/postprocess/hg38.py. It strips off those dotted
            # versions so that our resulting GTF file used by the workflows
            # will not contain them:

            postprocess: 'lib.postprocess.hg38.strip_ensembl_version'

            # Once well-formatted by the postprocessing function, we can now
            # perform standard conversions on the GTF. These conversions are
            # defined as rules in the references Snakefile, and will be run
            # if the conversion is specified here. Here we ask to get a refFlat
            # file, which can be provided to Picard's collectRnaSeqMetrics tool:

            conversions:
              - 'refflat'


        # Here is another tag, to create a FASTA file for ribosomal RNA. It can
        # then be used for fastq_screen, or for the rRNA screening portion of the
        # RNA-seq workflow:

        rRNA:
          fasta:

            # The SILVA database has separate files for large and small subunit
            # sequences. We'd like them all; by providing multiple URLs they will
            # be concatenated:

            url:
              - 'https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_LSURef_tax_silva_trunc.fasta.gz'
              - 'https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz'

            # However, the downloaded files contain many species. Here we only
            # care about human. We already have a function, "filter_fastas()", in
            # lib/common.py that accepts a FASTA and only keeps the records that
            # contain the provided first argument.

            # We specify that first argument here, and it will be passed to that
            # function, resulting in a final FASTA file that only contains the
            # rRNA sequence for Homo sapiens:

            postprocess:
                function: 'lib.common.filter_fastas'
                args: 'Homo sapiens'

            # We only need a bowtie2 index out of it.
            indexes:
                - 'bowtie2'

Without all those comments, it looks like this:

.. code-block:: yaml

    references_dir: 'data/references'
    organism: 'human'
    aligner:
      tag: 'gencode-v25'
      index: 'hisat2'
    references:
      human:
        gencode-v25:
          fasta:
            url: 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_organism.genome.fa.gz'
            indexes:
              - 'hisat2'
              - 'bowtie2'
          gtf:
            url: 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz'
            postprocess: 'lib.postprocess.hg38.strip_ensembl_version'
            conversions:
              - 'refflat'
        rRNA:
          fasta:
            url:
              - 'https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_LSURef_tax_silva_trunc.fasta.gz'
              - 'https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz'
            postprocess:
                function: 'lib.common.filter_fastas'
                args: 'Homo sapiens'
            indexes:
                - 'bowtie2'


The above file will result in the following directory structure::

    data/references/human/gencode-v25/fasta
    data/references/human/gencode-v25/bowtie2
    data/references/human/gencode-v25/hisat2
    data/references/human/gencode-v25/gtf
    data/references/human/gencode-v25-transcriptome/fasta
    data/references/human/gencode-v25-transcriptome/salmon
    data/references/human/rRNA/fasta
    data/references/human/rRNA/bowtie2

Each block in the YAML file describes either a `fasta` or `gtf` file. Each
block has at least the organism, type, and a URL.  A block can optionally have
a `postprocess`, which is an arbitrary function (described below) that converts
the downloaded URL to something that conforms to the standards of the workflow
(also described below). By supplying a tag, we can differentiate between
different versions (e.g., FlyBase r6.04 vs r6.11; hg19 vs hg38) or different
kinds of postprocessing (e.g, "chr" preprended to chrom names or not;
comprehensive annotation vs only coding genes).

`fasta` blocks can have an optional  `indexes` entry which will build the
specified indexes. `gtf` blocks can have an optional `conversions` entry which
will perform the specified conversion. Available indexes and conversions are
described below.


Post processing
---------------

**All files created by a block are required to be gzipped.**

This means that if a URL points to an uncompressed GTF file, a post-processing
function must gzip it. It also means that any post-processing functions must
write gzipped output files.

Other than that, it's up to the user to decide what transformations (if any)
are required. Examples might include:

- exluding particular contigs
- removing or editing problematic genes that have transcripts on both strands
  -- mod(mdg4) I'm looking at you
- renaming chromosomes (e.g., prepend "chr")
- remove unnecessary annotations (e.g., keep only cds/exon/transcript/gene features)

In the example config above, the yeast genome is available as a tarball of
separate fasta files, but we'd like to get it into a single fasta file for
downstream tools to work with.

The configuration block can define an optional `postprocess` string which
contains a dotted name referring to Python function that is importable by the
`reference.snakefile` workflow. By default, the workflow will find modules in
in ``lib.postprocess`` directory, so it's most convenient and organized to put
your functions within modules in that directory.

For example, above we used the postprocess function
``lib.postprocess.sacCer3.fasta_lib.postprocess``, and you can view this
function in ``lib/postprocess/sacCer3.py``.

Please see :func:`lib.common.download_and_postprocess` for more details, and
the files in the ``lib/postproces`` directory for inspiration.

These two arguments are automatically provided by the references workflow --
you don't have to know or care exactly what the filenames are, just what has to
be done to their contents.

See the files in ``lib/postprocess`` for inspiration if you need to write your
own post-processing functions.

The job of a postprocessing function is to ensure that the
fastq/gtf/transcriptome fasta meets the requirements described above and is
ready for any intended downstream tasks. For example if we download the fasta
file from FlyBase for dm6 but want "chr" prepended to chromosome names, we can
create a function in the file ``dm6.py`` called ``add_chr`` that does
this:

.. code-block:: python

    # This is dm6.py

    from snakemake.shell import shell  # a very convenient function

    def add_chr(origfn, newfn):
        shell(
            'zcat {origfn} '       # input is always gzipped
            '| sed "s/>/>chr/g" '  # add chr to names
            '| gzip -c > {newfn} ' # re-zip
            '&& rm {origfn}'       # clean up
        )

We specify this function to be called in the fasta config block like this (note
that the module doesn't have to be the same name as the organism, but it is
here for clarity):

.. code-block:: yaml

    dm6:
      fasta:
        url: ...
        postprocess: "dm6.add_chr"

This expects a file ``dm6.py`` in the same directory as the
`references.snakefile` workflow, and expects a function ``add_chr`` to
be defined in that module.

Any downstream rules that operate on the genome FASTA file (like hisat2 index,
bowtie2 index, etc) will now use this fixed version with "chr" prepended to
chromosome names.  In this way, we can apply arbitrary code to modify
references to get them into a uniform format.


More advanced postprocessing
----------------------------

If a post-processing function has a keyword argument with starts and ends with
a double underscore (``__``), the config system will assume this is a string
that should be interpreted as a dotted function name and the actual function
will be resolved and passed to the post-processing function.

This is useful for example when attaching ERCC spike-ins to a reference file
that in turn needs to be modified. For example, the `S. pombe` reference
annotations are available as a GFF file, but this needs to be converted to
a GTF file. After that, the ERCC spike-in GTF annotations need to be added to
the newly-created GTF.

The functions in ``lib/postprocess/ercc.py`` support such a use-case. The
config looks like this:

.. code-block:: yaml

    genome:
      url:
        # S. pombe fasta
        - 'ftp://ftp.ensemblgenomes.org/pub/fungi/release-41/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna_sm.toplevel.fa.gz'
        # ERCC fasta
        - 'https://www-s.nist.gov/srmors/certificates/documents/SRM2374_Sequence_v1.FASTA'
      postprocess:
        function: "lib.postprocess.ercc.add_fasta_to_genome"

    annotation:
      url:
        # S. pombe GFF, which needs to be converted to GTF
        - 'ftp://ftp.ensemblgenomes.org/pub/fungi/release-41/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.41.gff3.gz'

        # ERCC GTF is not available; conversion function needed to convert
        # fasta to GTF
        - 'https://www-s.nist.gov/srmors/certificates/documents/SRM2374_Sequence_v1.FASTA'

      postprocess:
        function: "lib.postprocess.ercc.add_gtf_to_genome"
        kwargs:
          # As per the docs for add_gtf_to_genome, this function will be
          # applied to all but the last input file. It is specified as a string
          # here, but the config-processing system will resolve this to the
          # actual function and pass that along to add_gtf_to_genome
          __preprocess__: "lib.common.gff2gtf"

.. versionadded:: 1.7
    Ability to use special ``__``-prefixed variables that are interpreted as
    dotted-path functions to import.

Locations of downloaded-and-post-rocessed FASTA and GTF files
-------------------------------------------------------------
Generally speaking, the fasta and gtf files will be in::

    {references_dir}/{organism}/{tag}/fasta/{organism}_{tag}.fasta
    {references_dir}/{organism}/{tag}/gtf/{organism}_{tag}.gtf

If a config file looks like this (simplified here for clarity):

.. code-block:: yaml

  references_dir: refs
  references:
    human:
      hg38:
        fasta: ...
        gtf: ...

Then the following files will be created::

    refs/human/hg38/fasta/human_hg38.fasta
    refs/human/hg38/gtf/human_hg38.gtf


If you are running the references workflow directly, or it is included in
another workflow that requests a chromsizes file, the following will also be
created::

    refs/human/hg38/fasta/human_hg38.chromsizes

.. note::

  URLs are expected to be gzipped and any postprocessing functions are
  expected to output gzipped files. This is because it is most common for
  providers to offer gzipped reference files, and therefore minimizes the
  effort required to prepare fasta and gtf files.  However, not all downstream
  tools handle gzipped input. The references workflow therefore stores only the
  uncompressed versions. We consider the resulting configuration simplicity to
  be worth the additional space and time cost.


Available indexes and conversions
---------------------------------
The following indexes can be currently be specified for fasta files:

hisat2
^^^^^^

    .. code-block:: yaml

        indexes:
          - hisat2

    Output files::

      {references_dir}/{organism}/{tag}/hisat2/{organism}_{tag}.*.ht2

bowtie2
^^^^^^^

    .. code-block:: yaml

        indexes:
          - bowtie2

    Output files::

      {references_dir}/{organism}/{tag}/bowtie2/{organism}_{tag}.*.bt2

salmon
^^^^^^

    .. code-block:: yaml

        indexes:
          - salmon

    Output files::

      {references_dir}/{organism}/{tag}/salmon/{organism}_{tag}/*

The following conversions can be specified for GTF files:

refflat
^^^^^^^

    .. code-block:: yaml

        conversions:
          - refflat

    Converts GTF to refFlat format. See the ``conversion_refflat`` rule in
    ``workflows/references/Snakefile``.

    Output file::

      {references_dir}/{organism}/{tag}/gtf/{organism}_{tag}.refflat

bed12
^^^^^

    .. code-block:: yaml

        conversions:
           - bed12

   Converts GTF to BED12 format. See the ``conversion_bed12`` rule in
   ``workflows/references/Snakefile``.

   Output file::

      {references_dir}/{organism}/{tag}/gtf/{organism}_{tag}.refflat

gffutils
^^^^^^^^
    Converts GTF to gffutils database (typically used for downstream work). You
    can specify arbitrary kwargs to ``gffutils.create_db`` by including them as
    keys. For example, if the GTF file already contains features for genes and
    transcripts:

    .. code-block:: yaml

        conversions:
          - gffutils:
              disable_infer_genes: True
              disable_infer_transcripts: True


    Output file::

        {references_dir}/{organism}/{tag}/gtf/{organism}_{tag}.gtf.db

genelist
^^^^^^^^
    Reads the postprocessed GTF file, and extracts the set of gene IDs found,
    one ID per line. The GTF attribute to use is configured by the
    ``gene_id:`` key, for example, if the file contains gene IDs in the
    ``Name`` attribute of each line, use the following:

    .. code-block:: yaml

        conversions:
          - genelist:
              gene_id: 'Name'

    Output file::

      {references_dir}/{organism}/{tag}/gtf/{organism}_{tag}.genelist

mappings
^^^^^^^^
    Reads the postprocesses GTF file, and outputs mappings between attributes
    as a gzipped TSV.

    You can include/exclude featuretypes from being checked.  For example, if
    your GTF has genes and transcripts in addition to exons, the gene and
    transcript lines probably contain all of the attributes you are interested
    in (like gene_id, symbol, name, etc) and the exon (and any other lines) can
    be ignored, speeding up the process. In this case you could use
    ``include_featuretypes: [gene, transcript]``.

    A ``__featuretype__`` column is always included in the mapping.  This is
    the GTF featuretype of each line, with extra ``__`` to avoid overwriting an
    attribute that may happen to be called ``featuretype``.

    .. code-block:: yaml

        conversions:
          - mappings

    .. code-block:: yaml

        conversions:
          - mappings:
              include_featuretypes: [gene, transcript]

    Output file::

      {references_dir}/{organism}/{tag}/gtf/{organism}_{tag}.mapping.tsv.gz
