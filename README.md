A streamlined version of LCDB workflows.

Guiding principles:

- Most of the work is in wrappers. The `lcdb-wrapper-tests` repo is included as
a submodule here.

- In the previous `lcdb-workflows`, the idea was to have everything in a single
config file. While great for end-users, it was too much overhead to write new
functionality in the workflows and then add all the other infrastructure in
order to expose that new functionality to the config file. Here, we allow
configuration to happen within the Snakefile (mostly via `params` fields) and
any config files remain lightweight.

- This package/set of workflows should also remain lightweight. Anything used
here should be familiar to anyone with Snakemake experience. There shouldn't be
any fancy infrastructure. Complexity should live in wrappers which should in
turn expose relatively simple APIs.

- The references workflow should ideally be run once per site; other workflows can
either point directly to the created files or can `include:` the workflow to
trigger updates

- Make heavy use of sampletables

- Any generally-useful helper functions go in `lcdblib`. Come to think of it we
may want to include that as submodule, too.

- It is expected that a particular workflow will get substantially edited
before actual use. Operating under the assumption that it's easier to delete
than to create, each workflow will have "the works" and can be trimmed down
according to the particular experiment's needs.

- Workflows should have a `patterns` dict at the top that lays out, in one
place, the output file patterns. Rules can optionally use values from this dict
to define output patterns. Using the `fill_patterns` function, these patterns
are "rendered" into a `targets` dictionary (basically a recursive `expand()`).
Selected contents of the `targets` directory are then used for the `all` rule.
This provides, all in one section of the Snakefile, a fair amount of
customization options. Patterns can be commented out, or targets can be
excluded from the `all` rule to fine-tune which rules will be run.

- A side effect of the `patterns` and `targets` design is that aggregation
rules become much easier to write. If all FastQC runs are under a `fastqc` key
in `targets`, they can all be used as input to a rule using
`flatten(targets['fastqc'])`. This is in contrast to, say, writing input
functions.


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
      - 'refflat'

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

Each block describes either a `fasta` or `gtf` file. Each block has at least
the assembly, type, and a URL.  They can also have a `postprocess`, which is an
arbitrary function (described below) that converts the downloaded URL to
something that conforms to the standards of the workflow (also described
below). By supplying a tag, we can differentiate between different versions
(e.g., FlyBase r6.04 vs r6.11) or different kinds of postprocessing (e.g, "chr"
preprended to chrom names or not).

Blocks with a type of `fasta` can have an optional  `indexes` entry which will
build the specified indexes. Blocks with a type of `gtf` can have an optional
`conversions` entry which will perform the specified conversion. Available
indexes and conversions are described below.

## Post processing
Sometimes the file at a URL is not in exactly the right format. The block can
define a `postprocess` string, which contains a dotted name referring to Python
function that is importable by the `reference.snakefile`. For example, if
a `postprocess` string is `"hg19.fix_fasta"`, then there should be a file
`hg19.py` that has within it a function called `fix_fasta()` in the same
directory as the references Snakefile. The dotted name should refer to
a function that has this function signature:

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

**All files created by a block are required to be gzipped.**

If a URL points to an uncompressed GTF file, a post-processing function must
gzip it. Any post-processing functions must write gzipped output files.

Other than that, it's up to the user to decide what transformations (if any)
are required. Examples might include:

* exluding contigs
* removing or editing problematic genes that have transcripts on both strands
* renaming chromosomes (e.g., prepend "chr")
* remove unnecessary annotations (e.g., keep only cds/exon/transcript/gene features)

## Available indexes

Currently:

    - hisat2
    - bowtie2
    - kallisto

Planned:

    - salmon
    - STAR
    - bwa

## Available conversions

Currently:

    - refflat (converts GTF to refFlat format)

Planned:

    - intergenic (needs chromsizes; therefor need to link a GTF tag to a FASTA
    tag but not quite sure how best to do this)
