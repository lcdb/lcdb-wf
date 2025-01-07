
rule multiqc:
    input:
        files=(
            expand(patterns["fastqc"]["raw"], sample=SAMPLES),
            expand(patterns["fastqc"]["cutadapt"], sample=SAMPLES),
            expand(patterns["fastqc"]["bam"], sample=SAMPLES),
            expand(patterns["markduplicates"]["bam"], sample=SAMPLES),
            expand(patterns["salmon"], sample=SAMPLES),
            expand(patterns["kallisto"], sample=SAMPLES),
            expand(patterns["preseq"], sample=SAMPLES),
            expand(patterns["rseqc"]["infer_experiment"], sample=SAMPLES),
            expand(patterns["rseqc"]["read_distribution"], sample=SAMPLES),
            expand(patterns["collectrnaseqmetrics"]["metrics"], sample=SAMPLES),
            expand(patterns["samtools"]["idxstats"], sample=SAMPLES),
            expand(patterns["samtools"]["flagstat"], sample=SAMPLES),
            expand(patterns["samtools"]["stats"], sample=SAMPLES),
            patterns["rrna_percentages_table"],
            patterns["featurecounts"],
        ),
        config='config/multiqc_config.yaml'
    output:
        'data/rnaseq_aggregation/multiqc.html'
    log:
        'data/rnaseq_aggregation/multiqc.log'
    threads: 1
    resources:
        mem_mb=gb(2),
        runtime=autobump(hours=2)
    run:
        analysis_directory = set([os.path.dirname(i) for i in input])
        outdir = os.path.dirname(output[0])
        basename = os.path.basename(output[0])
        shell(
            'LC_ALL=en_US.utf8 LC_LANG=en_US.utf8 '
            'multiqc '
            '--quiet '
            '--outdir {outdir} '
            '--force '
            '--filename {basename} '
            '--config {input.config} '
            '{analysis_directory} '
            '&> {log} '
        )
