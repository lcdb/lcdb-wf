Decision log
============

This document keeps track of the reasoning behind various architecture decisions.

References
----------
Here are use-cases we have that are common enough to warrant supporting:

- References should support multiple workflows (ChIP-seq, RNA-seq, etc)
  - This implies that the means the references dir should be in the
    ``workflows`` directory or above.
  - For example, this may mean a STAR index for RNA-seq, a bowtie2 index for
    rRNA contamination, and another bowtie2 index for ChIP-seq.

- References should support different organisms in different workflows. There
  should beo only one organism per workflow though.

- References should be re-created for each project.
  - What we've found is that if we have a central location for the references
    (shared by multiple deployments of lcdb-wf over the years) then we get
    conflicts where one deployment's aligner version is more recent, causing
    errors when using the index for an older version.
  - To keep using this, we'd need to version indexes based on aligner version.
  - However, when writing up methods for a paper we need to be able to trace
    back what commands were run to generate the reference, including additional
    patching that may have taken place (as is supported by the references
    workflow).
  - Re-using indexes is space- and time-efficient in the short term, but has
    shown to be inefficient in time and reproducibility in the long term.
  - Keeping everything in the same deployment director also helps with the
    archiving process. 

Naming:

- Top level should be organsim. Doesn't really matter in the case of
  a single-organism workflow.
- Next should be what has historically been called "tag". This could be the
  assembly name for genomic indexes, or some combination of assembly
  + annotation for transcriptome.
- If we're assuming "deployment-local" references, these no longer have to be
  globally unique. If we have a mouse reference with a transgene, we can just
  call it "mouse/mm39" but have the transgene patched into it, and not worry
  about conflicting (or worse, overwriting!) a central reference with the same
  name that didn't have the transgene.
- Fasta files are included next to their respective index.

This example uses the ``dmel`` organism and ``test`` tag which is configured by
default for tests.

This uses ``$ORG/$TAG/<genome|annotation|transcriptome>/$TOOL`` as the path
template. This lets us keep the fastq file used for building the various
indexes alongside the indexes.

::

  references_data/
  ├── dmel
      ├── rRNA
      │   └── genome
      │       ├── bowtie2
      │       │   └── dmel_rRNA.* <bowtie2 files>
      │       └── dmel_rRNA.fasta
      └── test
          ├── annotation
          │   ├── dmel_test.bed12
          │   ├── dmel_test.gtf
          │   └── dmel_test.refflat
          ├── genome
          │   ├── bowtie2
          │   │   └── dmel_test.* <bowtie2 files>
          │   ├── star
          │   │   └── dmel_test
          │   │       └── <STAR files>
          │   ├── dmel_test.chromsizes
          │   ├── dmel_test.fasta
          │   ├── dmel_test.fasta.fai
          └── transcriptome
              ├── kallisto
              │   └── dmel_test
              │       └── transcripts.idx
              ├── salmon
              │   └── dmel_test
              │       └── <salmon files>
              └── dmel_test.fasta
