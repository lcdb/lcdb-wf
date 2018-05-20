Config YAML
===========

============================================ =================== ================ ================= =========
Field                                        Used for References Used for RNA-seq Used for ChIP-seq Required
============================================ =================== ================ ================= =========
:ref:`references <cfg-references>`                    yes                 yes              yes      always
:ref:`references_dir <cfg-references-dir>`            yes                 yes              yes      if REFERENCES_DIR env var not set
:ref:`sampletable <cfg-sampletable>`                  .                   yes              yes      always
:ref:`organism <cfg-organism>`                        .                   yes              yes      always
:ref:`aligner <cfg-aligner>`                          .                   yes              yes      always
:ref:`gtf <cfg-gtf>`                                  .                   yes              .        always for RNA-seq
:ref:`rrna <cfg-rrna>`                                .                   yes              .        if rRNA screening desired
:ref:`salmon <cfg-salmon>`                            .                   yes              .        if Salmon quantification will be run
:ref:`fastq_screen <cfg-fastq-screen>`                .                   yes              yes      if using Fastq_screen
:ref:`merged_bigwigs <cfg-merged-bigwigs>`            .                   yes              yes      if you want to merge bigwigs
:ref:`chipseq <cfg-chipseq>`                          .                   .                yes      always for ChIP-seq
============================================ =================== ================ ================= =========

Example configs
---------------

RNA-seq
~~~~~~~

.. code-block:: yaml

    references_dir: "/data/references"
    sampletable: "config/sampletable.tsv"
    organism: 'human'
    aligner:
      tag: 'gencode-v25'
      index: 'hisat2'
    rrna:
      tag: 'rRNA'
      index: 'bowtie2'
    gtf:
      tag: '

    fastq_screen:
      - label: Human
        organism: human
        tag: gencode-v25
      - label: rRNA
        organism: human
        tag: rRNA

    # Portions have been omitted from "references" section below for
    # simplicity; see references config section for details.

    references:
      human:
        gencode-v25:
          fasta:
            url: 'ftp://.../genome.fa.gz'
            indexes:
              - 'hisat2'
              - 'bowtie2'
          gtf:
            url: 'ftp://.../annotation.gtf.gz'

        gencode-v25-transcriptome:
          fasta:
            url: 'ftp://.../transcriptome.fa.gz'
            indexes:
              - 'salmon'

        rRNA:
          fasta:
            url: 'https://...'
            indexes:
                - 'bowtie2'


Field descriptions
------------------
Required for references, RNA-seq and ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _cfg-references:

``references``
``````````````
    This section defines labels for references, where to get FASTA and GTF
    files and (optionally) post-process them, and which indexes to build. This
    is the most complex section; see :ref:`references-config` for details.

    Briefly, the example above has a single organism configured: "human". That
    organism has three tags: "gencode-v25", "gencode-v25-transcriptome", and
    "rRNA".

.. _cfg-references-dir:

``references_dir``
``````````````````
    Top-level directory in which to create references. If not specified the
    workflows will look for the environment variable ``REFERENCES_DIR``. If
    ``REFERENCES_DIR`` env var exists, it takes precedence over the
    ``references_dir`` field in the config file.

Required for RNA-seq and ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _cfg-sampletable:

``sampletable``
```````````````
    Path to sampletable file which, at minimum, list sample names and paths to
    FASTQ files. It is relative to the Snakefile. See :ref:`sampletable` for
    more info on the expected contents of the file.

    Example:

    .. code-block:: yaml

        sampletable: "config/sampletable.tsv"

.. _cfg-organism:

``organism``
````````````
    This field selects the top-level section of the ``references`` section that
    will be used for the analysis. In the example above, "human" is the only
    organism configured.

    Example:

    .. code-block:: yaml

        organism: "human"

.. _cfg-aligner:

``aligner``
```````````
    This field has two sub-fields, and automatically uses the configured
    ``organism`` to select the top-level entry in the references section.
    ``tag`` selects the tag from the organism to use, and ``index`` selects
    which aligner index to use. The relevant option from the example above
    would be "gencode-v25", which configures both bowtie2 and hisat2 indexes to
    be built. For RNA-seq we would likely choose "hisat2"; for ChIP-seq
    "bowtie2".

    Example:

    .. code-block:: yaml

        aligner:
          tag: "gencode-v25"
          index: "hisat2"

Optional fields
~~~~~~~~~~~~~~~

.. _cfg-fastq-screen:

``fastq_screen``
````````````````

    This section configures which Bowtie2 indexes should be used with
    `fastq_screen`. It takes the form of a list of dictionaries. Each
    dictionary has the keys:

        - `label`: how to label the genome in the output
        - `organism`: a configured organism. In the example above, there is only a single configured organism, "human".
        - `tag`: a configured tag for that organism.

    Each entry in the list must have a Bowtie2 index configured to be built.

    Example:

    .. code-block:: yaml

        fastq_screen:
          - label: Human
            organism: human
            tag: gencode-v25
          - label: rRNA
            organism: human
            tag: rRNA

.. _cfg-merged-bigwigs:

``merged_bigwigs``
``````````````````
    This section controls optional merging of signal files in bigWig format.
    Its format differs depending on RNA-seq or ChIP-seq, due to how strands are
    handled in those workflows.

    Here is an RNA-seq example:

    .. code-block:: yaml

        merged_bigwigs:
          arbitrary_label_to_use:
            sense:
              - 'sample1'
              - 'sample2'
            antisense:
              - 'sample1'
              - 'sample2'

    This will result in a single bigWig file called
    `arbitrary_label_to_use.bigwig` in the directory
    `data/rnaseq_aggregation/merged_bigwigs` (by default; this is configured
    using ``config/rnaseq_patterns.yaml``). That file merges together both the
    sense and antisense signal strands of two samples, sample1 and sample2. The
    names "sample1" and "sample2" are sample names defined in the :ref:`sample
    table <sampletable>`.

    Here's another RNA-seq example, where we merge the samples again but keep
    the strands separate. This will result in two output bigwigs.

    .. code-block:: yaml

        merged_bigwigs:
          merged_sense:
            sense:
              - 'sample1'
              - 'sample2'
          merged_antisense:
            antisense:
              - 'sample1'
              - 'sample

    Here is a ChIP-seq example:

    .. code-block:: yaml

        merged_bigwigs:
          arbitrary_label_to_use:
            - 'label1'
            - 'label2'

    This will result in a single bigWig file called
    `arbitrary_label_to_use.bigwig` in the directory
    `data/chipseq_aggregation/merged_bigwigs` (by default; this is configured
    using ``config/chipseq_patterns.yaml``) that merges together the "label1"
    and "label2" bigwigs.

    See :ref:`sampletable` for more info on the relationship between a *sample*
    and a *label* when working with ChIP-seq.


RNA-seq-only fields
~~~~~~~~~~~~~~~~~~~
.. _cfg-rrna:

``rrrna``
`````````

    This field selects the reference tag to use for screening rRNA reads.
    Similar to the ``aligner`` field, it takes both a ``tag`` and ``index``
    key. The specified index must have been configured to be built for the
    specified tag. It uses the already configured ``organism``.

    Example:

    .. code-block:: yaml

        rrna:
          tag: 'rRNA'
          index: 'bowtie2'


.. _cfg-gtf:

``gtf``
```````

    This field selects the reference tag to use for counting reads in features.
    The tag must have had a ``gtf:`` section specified; see
    :ref:`references-config` for details.

.. _cfg-salmon:

``salmon``
``````````
    This field selects the reference tag to use for the Salmon index (if used).
    The tag must have had a FASTA configured, and an index for "salmon" must
    have been configured to be built.

ChIP-seq-only fields
~~~~~~~~~~~~~~~~~~~~
.. _cfg-chipseq:

``chipseq``
```````````
    This section configures the peak-calling stage of the ChIP-seq workflow.
    This can get fairly complicated. 
