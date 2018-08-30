.. _config-yaml:

Config YAML
===========

Config files are expected to be in the file ``config/config.yaml`` relative to
the Snakefile. For example, the RNA-seq workflow at
``workflows/rnaseq/Snakefile`` expects the config file
``workflows/rnaseq/config/config.yaml``.

While it is possible to use Snakemake mechanisms such as ``--config`` to
override a particular config value and ``--configfile`` to update the config
with a different file, it is easiest to edit the existing
``config/config.yaml`` in place. This has the additional benefit of storing all
config information in one place for reproducibility.

The following table summarizes the config fields, which ones are use for which
workflow, and under what conditions, if any, they are required. Each option
links to a section below with more details on how to use it.

================================================================================== =================== ================ ================= =========
Field                                                                              Used for References Used for RNA-seq Used for ChIP-seq Required
================================================================================== =================== ================ ================= =========
:ref:`references <cfg-references>` and/or :ref:`include_references <cfg-inc-refs>`          yes                 yes              yes      yes
:ref:`references_dir <cfg-references-dir>`                                                  yes                 yes              yes      if `REFERENCES_DIR` env var not set
:ref:`sampletable <cfg-sampletable>`                                                        .                   yes              yes      always
:ref:`organism <cfg-organism>`                                                              .                   yes              yes      always
:ref:`aligner <cfg-aligner>`                                                                .                   yes              yes      always
:ref:`fastq_screen <cfg-fastq-screen>`                                                      .                   yes              yes      if using `fastq_screen`
:ref:`merged_bigwigs <cfg-merged-bigwigs>`                                                  .                   yes              yes      if you want to merge bigwigs
:ref:`gtf <cfg-gtf>`                                                                        .                   yes              .        always for RNA-seq
:ref:`rrna <cfg-rrna>`                                                                      .                   yes              .        if rRNA screening desired
:ref:`salmon <cfg-salmon>`                                                                  .                   yes              .        if Salmon quantification will be run
:ref:`chipseq <cfg-chipseq>`                                                                .                   .                yes      always for ChIP-seq
================================================================================== =================== ================ ================= =========

Example configs
---------------

To provide an overview, here are example config files.

RNA-seq
~~~~~~~

The config file for RNA-seq is expected to be in
``workflows/rnaseq/config/config.yaml``:

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
      tag: 'gencode-v25'

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

ChIP-seq
~~~~~~~~

The config file for ChIP-seq is expected to be in
``workflows/chipseq/config/config.yaml``.

The major differences between ChIP-seq and RNA-seq configs are:

- ChIP-seq has no ``gtf`` or ``rrna`` fields
- ChIP-seq has a ``chipseq: peak_calling:`` section

.. code-block:: yaml

    sampletable: 'config/sampletable.tsv'
    organism: 'dmel'

    aligner:
      index: 'bowtie2'
      tag: 'test'

    chipseq:
      peak_calling:

        - label: gaf-embryo-1
          algorithm: macs2
          ip:
            - gaf-embryo-1
          control:
            - input-embryo-1

        - label: gaf-embryo-1
          algorithm: spp
          ip:
            - gaf-embryo-1
          control:
            - input-embryo-1

        - label: gaf-wingdisc-pooled
          algorithm: macs2
          ip:
            - gaf-wingdisc-1
            - gaf-wingdisc-2
          control:
            - input-wingdisc-1
            - input-wingdisc-2

        - label: gaf-wingdisc-pooled
          algorithm: spp
          ip:
            - gaf-wingdisc-1
            - gaf-wingdisc-2
          control:
            - input-wingdisc-1
            - input-wingdisc-2

    fastq_screen:
      - label: PhiX
        organism: phix
        tag: default
      - label: Human
        organism: human
        tag: gencode-v25

    merged_bigwigs:
      input-wingdisc:
        - input-wingdisc-1
        - input-wingdisc-2
      gaf-wingdisc:
        - gaf-wingdisc-1
        - gaf-wingdisc-2
      gaf-embryo:
        - gaf-embryo-1


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

      dmel:
        test:
          fasta:
            url: "https://raw.githubusercontent.com/lcdb/lcdb-test-data/master/data/seq/dm6.small.fa"
            postprocess: 'lib.common.gzipped'
            indexes:
              - 'bowtie2'
              - 'hisat2'

      phix:
        default:
          fasta:
            url: 'ftp://...'
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
    files and (optionally) post-process them, and which indexes to build.

    Briefly, the example above has a single organism configured ("human"). That
    organism has three tags ("gencode-v25", "gencode-v25-transcriptome", and
    "rRNA").

    This is the most complex section and is documented elsewhere (see
    :ref:`references-config`).


.. _cfg-inc-refs:

``include_references``
``````````````````````

    This section can be used to supplement the ``references`` section with
    other reference sections stored elsewhere in files. It's a convenient way
    of managing a large amount of references without cluttering the config
    file.

    See :ref:`references-config` for more.


.. _cfg-references-dir:

``references_dir``
``````````````````
    Top-level directory in which to create references.

    If not specified, uses the environment variable ``REFERENCES_DIR``.

    If specified and ``REFERENCES_DIR`` also exists, ``REFERENCES_DIR`` takes
    precedence.

Required for RNA-seq and ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _cfg-sampletable:

``sampletable`` field
`````````````````````
    Path to sampletable file which, at minimum, list sample names and paths to
    FASTQ files. The path is relative to the Snakefile. See :ref:`sampletable`
    for more info on the expected contents of the file.

    Example:

    .. code-block:: yaml

        sampletable: "config/sampletable.tsv"

.. _cfg-organism:

``organism`` field
``````````````````
    This field selects the top-level section of the ``references`` section that
    will be used for the analysis. In the example above, "human" is the only
    organism configured.

    Example:

    .. code-block:: yaml

        organism: "human"

.. _cfg-aligner:

``aligner`` config section
``````````````````````````
    This field has two sub-fields, and automatically uses the configured
    ``organism`` to select the top-level entry in the references section.
    ``tag`` selects the tag from the organism to use, and ``index`` selects
    which aligner index to use. The relevant option from the example above
    would be "gencode-v25", which configures both bowtie2 and hisat2 indexes to
    be built. For RNA-seq we would likely choose "hisat2"; for ChIP-seq
    "bowtie2".

    Currently-configured options are ``hisat2``, ``bowtie2``, ``star``, and
    ``ngm``.

    Example:

    .. code-block:: yaml

        aligner:
          tag: "gencode-v25"
          index: "hisat2"

Optional fields
~~~~~~~~~~~~~~~

.. _cfg-fastq-screen:

``fastq_screen`` config section
```````````````````````````````

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

``merged_bigwigs`` config section
`````````````````````````````````
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

``rrrna`` field
```````````````

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

``gtf`` field
`````````````

    This field selects the reference tag to use for counting reads in features.
    The tag must have had a ``gtf:`` section specified; see
    :ref:`references-config` for details.

.. _cfg-salmon:

``salmon`` field
````````````````
    This field selects the reference tag to use for the Salmon index (if used).
    The tag must have had a FASTA configured, and an index for "salmon" must
    have been configured to be built for the organism selected with the
    ``organism`` config option.

ChIP-seq-only fields
~~~~~~~~~~~~~~~~~~~~

.. _cfg-chipseq:

``chipseq`` config section
``````````````````````````
    This section configures the peak-calling stage of the ChIP-seq workflow. It
    currently expects a single key, ``peak_calling``, which is a list of
    peak-calling runs.

    A peak-calling run is a dictionary configuring a single peak-calling run
    which results in a single BED file of called peaks. A peak-calling run is
    uniquely described by its ``label`` and ``algorithm``. This way, we can use
    the same label (e.g., `gaf-embryo-1`) across multiple peak-callers to help
    organize the output.

    Here is a minimal example of a peak-calling config section. It defines
    a single peak-calling run using the `macs2` algorithm. Note that the
    ``ip:`` and ``control:`` keys are lists of **labels** from the ChIP-seq
    sample table's ``label`` column, not sample IDs from the first column.

    .. code-block:: yaml

        chipseq:
          peak_calling:

            - label: gaf-embryo-1
              algorithm: macs2
              ip:
                - gaf-embryo-1
              control:
                - input-embryo-1

    The above peak-calling config will result in a file
    ``data/chipseq_peaks/macs2/gaf-embryo-1/peaks.bed`` (that pattern is
    defined in ``chipseq_patterns.yaml`` if you need to change it).

    We can specify additional command-line arguments that are passed verbatim
    to `macs2` with the ``extra:`` section, for example:

    .. code-block:: yaml

        chipseq:
          peak_calling:

            - label: gaf-embryo-1
              algorithm: macs2
              ip:
                - gaf-embryo-1
              control:
                - input-embryo-1
              extra: '--nomodel --extsize 147'


    `macs2` supports multiple IP and input files, which internally are merged
    by `macs2`. We can supply multiple IP and input labels for biological
    replicates to get a set of peaks called on pooled samples. Note that we
    give it a different label so it doesn't overwrite the other peak-calling
    run we already have configured.

    .. code-block:: yaml

        chipseq:
          peak_calling:

            - label: gaf-embryo-1
              algorithm: macs2
              ip:
                - gaf-embryo-1
              control:
                - input-embryo-1
              extra: '--nomodel --extsize 147'


            - label: gaf-embryo-pooled
              algorithm: macs2
              ip:
                - gaf-embryo-1
                - gaf-embryo-2
              control:
                - input-embryo-1
                - input-embryo-2


