
.. _references-config:

References config
=================

The references section defines genomes, transcriptomes, and annotations to use.
It supports arbitrarily many species and assemblies, and supports customizing
references for a particular project.

Using existing configs
----------------------
If you just want to use some references configs that work, then put this in
your config file:

.. code-block:: yaml

    include_references:
      - '../../include/references_configs'

This will populate the config with the contents of all the files contained in
``include/references_configs``; the path configured here is relative to the
Snakefile using the config. You will still need to inspect the contents of
those files to decide which organsim and tag you want to use for your
particular experiment (see :ref:`cfg-organism` and :ref:`cfg-aligner` for more
on these fields). For example, if you are working with human RNA-seq data, you
may want this in your config:

.. code-block:: yaml

    organism: 'human'
    aligner:
      tag: 'gencode-v25'
      index: 'hisat2'
    salmon: 'gencode-v25-transcriptome'

The remainder of this section of the documentation explains how to customize
the references, to add your own or modify the existing examples.

Overview
--------
The references workflow is based on the idea that while each genome's source
files (FASTA, GTF) may come from different providers (Ensembl, FlyBase, UCSC,
etc) and have slightly different formatting (strange headers, one big file or
a tarball of individual chromosomes, etc), once they are well-formatted they
can be used to create a hisat2 index, a bowtie2 index, a list of genes,
intergenic regions, and so on without any further customization.

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

It's probably easiest to show an example config and then describe what's
happening.

The following example configures the workflow to download a tarball of fasta
files for the yeast genome from UCSC, concatenate them into a single file, and
build a bowtie2 and hisat2 index for it. It's heavily commented for
illustration.


.. code-block:: yaml

    # This configures the directory in which the prepared references will be
    # saved (see below for directory structure):

    references_dir: 'data/references'

    # One of the organisms configured below. We are only configuring a single one
    # so sacCer3 is our only option here:

    organism: 'saccharomyces_cerevisiae'

    # Here we specify which tag (under organism "sacCer3" below) to use for
    # aligning, as well as which index we'll be using. This example is RNA-seq,
    # so we'll use HISAT2:

    aligner:
      tag: 'sacCer3'
      index: 'hisat2'

    # Top-level section for references:

    references:

      # Label for this organism or species:

      saccharomyces_cerevisiae:

        # "sacCer3" is our tag to describe this particular FASTA and GTF we're
        # preparing:

        sacCer3:

          # This block will define how to get and postprocess a FASTA file:

          fasta:

            # URL to download:

            url: 'http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz'

            # The yeast genome from UCSC comes as a tarball of fastas. We can
            # specify a function to apply to the downloaded tarball to get
            # a single fasta file. See below for details.

            postprocess: 'lib.postprocess.sacCer3.fasta_lib.postprocess'

            # We can optionally build indexes for various aligners:

            indexes:
                - 'bowtie2'
                - 'hisat2'

Without all those comments, it looks like this:

.. code-block:: yaml

    references_dir: 'data/references'
    organism: 'saccharomyces_cerevisiae'
    aligner:
      tag: 'sacCer3'
      index: 'hisat2'
    references:
      saccharomyces_cerevisiae:
        sacCer3:
          fasta:
            url: 'http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz'
            postprocess: 'lib.postprocess.sacCer3.fasta_lib.postprocess'
            indexes:
                - 'bowtie2'
                - 'hisat2'

Each block in the YAML file describes either a `fasta` or `gtf` file. Each
block has at least the organism, type, and a URL.  A block can optionally have
a `postprocess`, which is an arbitrary function (described below) that converts
the downloaded URL to something that conforms to the standards of the workflow
(also described below). By supplying a tag, we can differentiate between
different versions (e.g., FlyBase r6.04 vs r6.11) or different kinds of
postprocessing (e.g, "chr" preprended to chrom names or not).

Blocks with a type of `fasta` can have an optional  `indexes` entry which will
build the specified indexes. Blocks with a type of `gtf` can have an optional
`conversions` entry which will perform the specified conversion. Available
indexes and conversions are described below.

Running the references workflow using this config will result in the following
files in ``data/references``::

      saccharomyces_cerevisiae
      └── sacCer3
          ├── bowtie2
          │   ├── saccharomyces_cerevisiae_sacCer3.1.bt2   # bowtie2 index files
          │   ├── saccharomyces_cerevisiae_sacCer3.2.bt2
          │   ├── saccharomyces_cerevisiae_sacCer3.3.bt2
          │   ├── saccharomyces_cerevisiae_sacCer3.4.bt2
          │   ├── saccharomyces_cerevisiae_sacCer3.rev.1.bt2
          │   └── saccharomyces_cerevisiae_sacCer3.rev.2.bt2
          └── fasta
              ├── saccharomyces_cerevisiae_sacCer3.chromsizes
              ├── saccharomyces_cerevisiae_sacCer3.fasta
              └── saccharomyces_cerevisiae_sacCer3.fasta.gz.log


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

A post-processing function has the following signature:

.. code-block:: python

    def func(temp_downloaded_filenames, final_postprocessed_filename, **kwargs)

where:

    - ``temp_downloaded_filenames`` is a list of temporary files downloaded for
      each provided URL. In our example above, there's only one URL provided,
      so this will be a list of one item.
    - ``final_postprocessed_filename`` is the final filename to create and will
      be a string.

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


TODO: document the conversions for GTF, specifically the `genelist` and
`annotation_hub` conversions and how the kwargs can be specified.

Available indexes and conversions
---------------------------------

Current indexes:

    - hisat2
    - bowtie2
    - kallisto
    - salmon

Planned indexes:

    - STAR
    - bwa

Current conversions:

    - refflat (converts GTF to refFlat format)
    - gffutils (converts GTF to gffutils database)
    - genelist
    - annotation_hub

Planned:

    - intergenic (needs chromsizes; therefore need to link a GTF tag to a FASTA
      tag but not quite sure how best to do this)
