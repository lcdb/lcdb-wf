.. _chipseq:

ChIP-seq workflow
-----------------
The ChIP-seq workflow includes:

    - trimming reads with cutadapt
    - mapping reads with Bowtie2
    - FastQC on raw, trimmed, and aligned reads
    - Remove multimappers (samtools) and duplicates (Picard MarkDuplicates)
    - fastq_screen on multiple configured genomes to look for evidence of
      cross-contamination
    - QC aggregation using MultiQC, along with a custom table for library sizes
    - merging and re-deduplicating to correctly handle technical replicates
    - bigWigs created from unique, no-dups BAM files
    - deepTools plotFingerprint run on grouped IP and input for QC and
      evaluation of enrichment
    - peak-calling using macs2 and/or spp
    - conversion of BED files into bigBed (or bigNarrowPeak where possible)
    - track hub of bigWigs and bigBeds to visualize peak-calling in UCSC Genome Browser

.. image:: chipseq.png
