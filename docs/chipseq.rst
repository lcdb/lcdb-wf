.. _chipseq:

ChIP-seq workflow
-----------------
The ChIP-seq workflow starts with raw FASTQ files and performs various QC steps. It
aligns and prepares BAM and bigWig files, performs peak-calling, and combines
everything together into a track hub for visualization.

Specifically, the workflow does the following:

    - trim reads with cutadapt
    - map reads with Bowtie2
    - run FastQC on raw, trimmed, and aligned reads
    - Remove multimappers (samtools) and duplicates (Picard MarkDuplicates)
    - fastq_screen on multiple configured genomes to look for evidence of
      cross-contamination
    - QC aggregation using MultiQC, along with a custom table for library sizes
    - merge technical replicates and then re-deduplicate
    - create bigWigs from unique, no-dups BAM files
    - optionally merge bigWigs to create one signal track for all replicates
    - run deepTools plotFingerprint on grouped IP and input for QC and
      evaluation of enrichment
    - call peaks using macs2, spp, and/or sicer, with support for multiple
      peak-calling runs using different parameters to assist with assessing
      performance and to help make decisions for downstream analysis
    - optionally run a template diffBind RMarkdown file used for differential binding analysis
    - convert BED files into bigBed (or bigNarrowPeak where possible)
    - build and optionally upload a track hub of bigWigs and bigBeds to
      visualize peak-calling in UCSC Genome Browser

To configure a ChIP-seq experiment, see :ref:`config-yaml`.

.. image:: chipseq.png
