.. _references:

References workflow
===================

This workflow is intended to be `include:`-ed into other workflows that depend
on reference fastas, indexes, and annotations. In that case, rules in this
references workflow will only be run for those files asked for in the parent
workflow.

This workflow can also be called on its own, in which case it will build
**all** of the references and indexes specified in the config. This can be
helpful when setting up the workflows for the first time on a new machine. Run
it like this::

    snakemake -s references.snakefile --configfile config/config.yaml

In both cases, it depends on the `references` section being in the global
config dictionary. See below for details.

Overview
--------
The references workflow is based on the idea that while each genome's source
files (FASTA, GTF) may come from different places and have slightly different
formatting, once they are well-formatted they can be used to create a hisat2
index, a list of genes, intergenic regions, and so on without any further
customization.

The challenging part is the "well-formatted" part. The config file allows very
flexible specification of how to create the files by providing a sort of plugin
architecture.

It's probably easiest to show an example config and then describe what's
happening. This example downloads a tarball of fasta files for the yeast
genome, concatenates them into a single and builds a bowtie2 index for it.
A chromsizes file is also built.


.. code-block:: yaml

    references_dir: 'data/references'

    # top-level section for references
    references:

      # label for the assembly or species
      sacCer3:

        # "tag" to differentiate between different versions.
        default:

          # this block will define how to get and postprocess a FASTA file
          fasta:

            # URL to download
            url: 'http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz'

            # The yeast genome from UCSC comes as a tarball of fastas. We can
            # specify a function to apply to the downloaded tarball to get
            # a single fasta file. See below for details.
            postprocess: 'lib.postprocess.sacCer3.fasta_lib.postprocess'

            # We can optionally build indexes for various aligners.
            indexes:
                - 'bowtie2'
                - 'hisat2'

Each block in the YAML file describes either a `fasta` or `gtf` file. Each
block has at least the assembly, type, and a URL.  A block can optionally have
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

      sacCer3/
      └── default
          ├── bowtie2
          │   ├── sacCer3_default.1.bt2   # bowtie2 index files
          │   ├── sacCer3_default.2.bt2
          │   ├── sacCer3_default.3.bt2
          │   ├── sacCer3_default.4.bt2
          │   ├── sacCer3_default.rev.1.bt2
          │   └── sacCer3_default.rev.2.bt2
          └── fasta
              ├── sacCer3_default.chromsizes
              ├── sacCer3_default.fasta
              └── sacCer3_default.fasta.gz.log


Requirements
------------

**All files created by a block are required to be gzipped.**

If a URL points to an uncompressed GTF file, a post-processing function must
gzip it. Any post-processing functions must write gzipped output files.

Other than that, it's up to the user to decide what transformations (if any)
are required. Examples might include:

* exluding particular contigs
* removing or editing problematic genes that have transcripts on both strands
  -- mod(mdg4) I'm looking at you
* renaming chromosomes (e.g., prepend "chr")
* remove unnecessary annotations (e.g., keep only cds/exon/transcript/gene features)

Post processing
---------------

Sometimes the file at a URL is not in exactly the right format. In the example
above, the yeast genome is available as a tarball of separate fasta files, but
we'd like to get it into a single fasta file for downstream tools to work with.

The configuration block can define an optional `postprocess` string which
contains a dotted name referring to Python function that is importable by the
`reference.snakefile` workflow.  For example, if a `postprocess` string is
`"hg19.fix_fasta"`, then there should be a file `hg19.py` that has within it
a function called `fix_fasta()` in the same directory as the references
Snakefile. The dotted name should refer to a function that has this function
signature:

.. code-block:: python

    def func(temp_downloaded_filenames, final_postprocessed_filename)


The first argument is a list corresponding to the tempfiles downloaded for each
provided url; the second is the final filename to create. These two arguments
are automatically provided by the references workflow -- you don't have to know
or care exactly what the filenames are, just what has to be done to their
contents.

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
that the module doesn't have to be the same name as the assembly, but it is
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
