lcdb-wf is a collection of snakemake workflows and tools for common high-throughput sequencing analysis, along with associated infrastructure.

See docs at https://lcdb.github.io/lcdb-wf.

RNASEQ workflow

This workflow is used for RNA-seq and RNA-seq-like analysis (like euRNA-seq, RIP-seq or small RNA-seq).

This workflow can use references created by the references workflow with no need to run the references workflow separately. This workflow performs the following tasks:

- Build a HISAT2 index
- Build a salmon transcriptome index
- Download a GTF annotation
- Convert the GTF to refflat format
- Trim reads with cutadapt
- Align with HISAT2
- Run FastQC on raw, trimmed, and aligned reads
- Align reads to rRNA using bowtie2 to evaluate rRNA contamination
- Count reads in genes with featureCounts
- Run dupRadar and preseq to assess library complexity
- Check for evidence of cross-contamination using fastq_screen on multiple configured genomes
- Assess transcript coverage with Picard CollectRnaSeqMetrics
- Build bigWigs (optionally strand-specific) created from BAM files
- Optionally merge bigWigs as defined by config
- Aggregate QC results using MultiQC. Includes custom tables for library sizes and rRNA contamination
- Run various QC and differential expression. This is performed in an RMarkdown file that runs a standard DESeq2 differential expression analysis along with diagnostic plots, exported tables of differentially expressed genes for each contrast, and downstream GO analysis using clusterProfiler. This file is run and rendered into an output HTML file.
- Construct and upload a track hub of scaled coverage bigWigs for each sample that can be viewed in UCSC Genome Browser


Configurations used:

{% for items in snakemake.config %}
  {% if items != 'references'%}
    - {{ items }} : {{snakemake.config[items]}}
  {% elif items == 'references' %}
    - References:

    {% for sublista in snakemake.config[items] %}
      - {{sublista}} :

      {% for sublistb in snakemake.config[items][sublista] %}
        - {{ sublistb }} :

        {% for sublistc in snakemake.config[items][sublista][sublistb] %}

          - {{sublistc}} : {{snakemake.config[items][sublista][sublistb][sublistc]}}

        {% endfor %}
      {% endfor %}
    {% endfor %}
  {% endif %}
{% endfor %}


