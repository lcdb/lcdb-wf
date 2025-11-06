.. _workflows:

Workflows
=========

The workflows are RNA-seq-like, ChIP-seq-like, and variant calling.

They are currently labeled ``rnaseq``, ``chipseq``, and ``variant-calling``,
but you can rename the workflow directories to whatever you need. And you can
make copies of a workflow directory to support multiple experiments.

If multiple experiments can all use the same parameters in the Snakefile, all
the samples can be combined into the same sampletable. But if they differ --
for example, they had different library prep such that the cutadapt parameters
need to be changed -- then they need to split up into multiple workflow
directories, each with their own sample and with respective edits to the
Snakefile. See :ref:`decisions-sample-specific-params` for rationale.

RNA-seq
-------

**Sampletable:** :ref:`rnaseq-sampletable`

**Config:** :ref:`rnaseq-config`

**Downstream:** :ref:`rnaseq-downstream`

This workflow can be used for any bulk Illumina-based RNA-seq-like assay that
quantifies transcripts of some sort and where a gene-by-sample matrix of counts
is useful.

This of course includes standard bulk RNA-seq, but also things like
RIP-seq, small RNA-seq, or even differential ChIP-seq within gene
bodies.

This workflow trims raw reads with cutadapt, aligns with STAR and quantifies
reads in genes with featureCounts. It also quantifies reads with Salmon.
Extensive QC is performed at each stage and is aggregated with MultiQC.

The biggest advantage of using this workflow is the extensive downstream
analysis (see :ref:`rnaseq-downstream`), which is run after the Snakefile
completes due to the frequent need for project-specific customization.

The primary output of the Snakefile consists of the following:

- Salmon quantification files for each sample

::

  data/rnaseq_samples/{sample_id}/{sample_id}.salmon/quant.sf


- aligned BAM files for each sample (duplicates marked but not removed)

::

   data/rnaseq_samples/{sample_id}/{sample_id}.cutadapt.markdups.bam

- strand-specific bigWig files for each sample

::

  data/rnaseq_samples/{sample_id}/{sample_id}.cutadapt.bam.neg.bigwig
  data/rnaseq_samples/{sample_id}/{sample_id}.cutadapt.bam.pos.bigwig

- single featureCounts file with all samples

::

  data/rnaseq_aggregation/featurecounts.txt

- MultiQC output

::

  data/rnaseq_aggregation/multiqc.html

The primary output of the downstream analysis (:ref:`rnaseq-downstream`) is the
final HTML report and the RDS files ready for exploration with `Carnation
<https://github.com/NICHD-BSPC/carnation>`__.

::

   downstream/rnaseq.html
   downstream/combined.Rds

ChIP-seq (and other chromatin-associated assays)
------------------------------------------------

**Sampletable:** :ref:`chipseq-sampletable`

**Config** :ref:`chipseq-config`

This workflow can be used for various bulk Illumina-based sequencing assays related
to chromatin binding, like ChIP-seq, CUT&RUN, Cut&Tag, ATAC-seq. There may need
to be modifications you need to make within the particular tool calls, but the
framework is useful for all of them.

This workflow trims raw reads with cutadapt, aligns with bowtie2, and runs peak
calling on all samples. Extensive QC is performed at each stage which is
aggregated with MultiQC.

The biggest advantage of using this workflow is the flexibility of
peak-calling. Since peak-calling tends to need extensive tweaking depending on
the antibody or assay, it is straightfoward to configure multiple peak-calling
runs (different algorithmss, each with possibley different parameters) on the
same sample, and view them all together in a genome browser to decide on
a final strategy.

The primary output of the Snakefile consists of the following:

- aligned BAM files (multimappers removed, duplicates removed)

::

  data/chipseq_samples/{sample_id}/{sample_id}.cutadapt.unique.nodups.bam

- peak calls

::

  data/chipseq_peaks/{algorithm}/{peak_run}/peaks.bed
  data/chipseq_peaks/{algorithm}/{peak_run}/peaks.bigbed

- bigWigs for merged technical replicates

::

  data/chipseq_merged/{biological_material}/{biological-material}.cutadapt.unique.nodups.bam.bigwig

- MultiQC output

::

  data/chipseq_aggregation/multiqc.html


Multiple workflows
------------------

If you have multiple experiments of the same type, and the same parameters can
be used for all of them, then they can go in the same sampletable. 

If samples need different parameters in any way, then make a copy of the
respective workflow directory and consider them as part of a different workflow.

(see :ref:`decisions-sample-specific-params` for rationale on this)

For example, miRNA-seq and SMART-seq would likely need different cutadapt
parameters and maybe alignment parameters, so we might put them in different
workflows:

::

  workflows/
    mirnaseq/
    smartseq/

Each of these would be a copy of the ``rnaseq`` workflow, but with appropriate
changes in the respective Snakefiles.

This is also the mechanism for working with different genomes, which would have different references in the config::

  workflows/
    chipseq-human/
    chipseq-mouse/

Each workflow can be considered as independent, which gives lots of flexibility
in configuring and customizing the Snakefile.
