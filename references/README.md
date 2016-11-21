# References workflow

This workflow can be called from other workflows that depend on reference
fastas, indexes, and annotations. In that case, rules in this workflow will
only be run for those files asked for in the parent workflow.

This workflow can also be called on its own, in which case it will build
**all** of the references and indexes specified in the config. This can be
helpful when setting up the workflows for the first time on a new machine. Run
it like this:

```bash
snakemake -s references.snakefile --configfile ../mapping/config.yaml
```

In both cases, it depends on the `references` section being in the global
config dictionary; see below for details.

## Overview

The references workflow is intended to be a universal workflow that supports
arbitrary genomes and results in FASTAs, indexes, annotations, and various
processed versions of these files. It is based on the idea that while each
genome's source files (FASTA, GTF) may come from different places and have
slightly different formatting, once they are well-formatted they can be used to
create a hisat2 index, a list of genes, intergenic regions, and so on without
any further customization.

The config file allows very flexible specification of how to
create the files by providing a sort of plugin architecture.

It's probably easiest to show an example config and then describe what's
happening. This example sets up a complete set of Drosophila references and
indexes, and just a bowtie2 index for E. coli (for use by fastq_screen to check
for contamination).

```yaml

# The following config is part of the main config.yaml, rather than
# a standalone file

data_dir: /data/LCDB/references
references:
  -
    assembly: 'dm6'
    type: 'gtf'
    url: 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gtf/dmel-all-r6.11.gtf.gz'
    postprocess: "dm6.gtf_postprocess"
    tag: 'r6-11'
    conversions:
      - 'intergenic'

  -
    assembly: 'dm6'
    url: 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.11.fasta.gz'
    postprocess: "dm6.fasta_postprocess"
    type: 'fasta'
    indexes:
      - 'bowtie2'
      - 'hisat2'

  -
    assembly: 'dm6'
    type: 'fasta'
    tag: 'r6-11_transcriptome'
    url: 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-transcript-r6.11.fasta.gz'
    indexes:
      - 'kallisto'

  -
    assembly: 'ecoli'
    type: 'fasta'
    url: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz'
    indexes:
      - 'bowtie2'

  -
    assembly: 'wolbachia'
    type: 'fasta'
    url: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.fna.gz'
    indexes:
      - 'bowtie2'

```

Each block describes either a fasta or gtf file. Each block has at least the
assembly, type, and a URL.  They can also have a `postprocess`, which is an
arbitrary function (described below) that converts the downloaded URL to
something that conforms to the standards of the workflow (also described
below). By supplying a tag, we can differentiate between different versions
(e.g., FlyBase r6.04 vs r6.11) or different kinds of postprocessing (e.g, "chr"
preprended to chrom names or not).

Blocks with a type of "fasta" can have an optional  `indexes` entry which will
build the specified indexes.

## Post processing
The `postprocessing` values are dotted names that refer to Python modules
importable by the `reference.snakefile`. The dotted name should refer to
a function that has the function signature:

```python
def func(temp_downloaded_filenames, final_postprocessed_filename)
```

These two arguments are automatically provided by the references workflow --
you don't have to know or care exactly what the filenames are, just what has to
be done to their contents. The first argument is a list corresponding to the
tempfiles downloaded for each provided url; the second is the final filename to
create.

The job of a postprocessing function is to ensure that the
fastq/gtf/transcriptome fasta meets the requirements below and is ready for any
intended downstream tasks. For example if we download the fasta file from
FlyBase for dm6 but want "chr" prepended to chromosome names, we can create
a function in the file `dm6.py` called `fasta_postprocess` that does this:

```python
def fasta_postprocess(origfn, newfn):
    shell("""zcat {origfn} | sed "s/>/>chr/g" | gzip -c > {newfn}  && rm {origfn}""")
```

Note that the file is uncompressed in order to fix the chromosome naming and
then recompressed into the final output. In this case the original file is
removed.  We specify this function to be called in the fasta config block like
this (note that the module doesn't have to be the same name as the assembly,
but it is here for clarity):

```yaml
dm6:
  fasta:
    url: ...
    postprocess: "dm6.fasta_postprocess"
```

This expects a file `dm6.py` in the same directory as the
`references.snakefile` workflow, and expects a function `fasta_postprocess` to
be defined in that module.

Any downstream rules that operate on the genome FASTA file (like hisat2 index,
bowtie2 index, etc) will now use this fixed version with "chr" prepended to
chromosome names.  In this way, we can apply arbitrary code to modify
references to get them into a uniform format.

## Requirements

All files are required to be gzipped. Note that there is a `unzip_fasta` rule
that will temporarily unzip the fasta for things like bowtie2 index building
that don't support gzipped input.

Other than that, it's up to the user to decide what transformations (if any)
are required. Examples might include:

* exluding contigs
* removing or editing problematic genes that have transcripts on both strands
* renaming chromosomes (e.g., prepend "chr")
* remove unnecessary annotations (e.g., keep only cds/exon/transcript/gene features)

## TODO

A rule should be triggered to be re-run based on md5sum of a concatenation of
the URL and the code object of the post-processing function.

How to differentiate between flybase r6.06 and r6.11? There should be a *tag*
associated with each config block. Maybe borrow the build string idea from
conda -- `{assembly}-{tag}.{extension}` or something.

How to handle cases like the hg19 GENCODE, where protein-coding transcripts are
in one FASTA but lncRNAs are in another?
