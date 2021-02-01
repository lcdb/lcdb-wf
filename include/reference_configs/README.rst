The files in this directory provide example configuration files for various
genomes (or sets of genomes).

Paste their contents into the ``references:`` section of a config file you're
using for the references, RNA-seq, or ChIP-seq workflow (typically in
``config/config.yaml``). Note that when doing this, do not paste the first
``references:`` line.

Alternatively, provide one or more config files to Snakemake with the
``--configfile`` argument. Note that this will *update* any existing contents
of the ``references:`` section of the default config file.

If you want to generate all references provided here, create one big file with
everything, like this:

.. code-block:: bash

    echo "references:" > all.yaml && cat *yaml | grep -v "^references:" >> all.yaml

This can be useful when setting up a new site for the first time.

For example, consider the following config file in ``config/config.yaml``:

.. code-block:: yaml

    sampletable: 'config/sampletable.tsv'
    organism: 'genomeA'
    references_dir: 'references_data'

    ... (other config items) ...

    references:
      genomeA:
        default:
          fasta:
            url: 'http://mygenome.org/a.fa.gz'
            indexes:
              - 'bowtie2'


and the following separate references file, in
``../../include/references_configs/GenomeB.yaml``:

.. code-block:: yaml

    references:
      genomeB:
        default:
          fasta:
            url: 'http://mygenome.org/b.fa.gz'
            indexes:
              - 'bowtie2'
        transcriptome:
          fasta:
            url: 'http://anothergenome.org/t.fa.gz'
            indexes:
              - 'salmon'

Then the following command would update the default ``config/config.yaml``'s
references section with genome B's information, allowing both genome A and
genome B to be available to the Snakefile:


.. code-block:: bash

    cd lcdb-wf/workflows/references
    snakemake --configfile ../../include/references_configs/GenomeB.yaml --use-conda -j8

Doing so effectively gives the Snakefile the following config file:

.. code-block:: yaml

    sampletable: 'config/sampletable.tsv'
    organism: 'genomeA'
    references_dir: 'references_data'

    ... (other config items) ...

    references:
      genomeA:
        default:
          fasta:
            url: 'http://mygenome.org/a.fa.gz'
            indexes:
              - 'bowtie2'

      # These are the contents from
      # ../../include/references_configs/GenomeB.yaml
      genomeB:
        default:
          fasta:
            url: 'http://mygenome.org/b.fa.gz'
            indexes:
              - 'bowtie2'
        transcriptome:
          fasta:
            url: 'http://anothergenome.org/t.fa.gz'
            indexes:
              - 'salmon'
