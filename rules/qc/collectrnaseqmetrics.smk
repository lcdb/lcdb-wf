rule collectrnaseqmetrics:
    """
    Calculate various RNA-seq QC metrics with Picarc CollectRnaSeqMetrics
    """
    input:
        bam=patterns['markduplicates']['bam'],
        refflat=rules.conversion_refflat.output,
    output:
        metrics=patterns['collectrnaseqmetrics']['metrics'],
    params:
        java_args='-Xmx20g',
        # java_args='-Xmx2g',  # [TEST SETTINGS -1]
        strand_arg={
                'unstranded': 'STRAND=NONE ',
                'fr-firststrand': 'STRAND=SECOND_READ_TRANSCRIPTION_STRAND ',
                'fr-secondstrand': 'STRAND=FIRST_READ_TRANSCRIPTION_STRAND ',
            }[config["stranded"]]
    log:
        patterns['collectrnaseqmetrics']['metrics'] + '.log'
    threads: 1
    resources:
        mem_mb=gb(32),
        runtime=autobump(hours=2)
    run:
        shell(
            'picard '
            '{params.java_args} '
            'CollectRnaSeqMetrics '
            '{params.strand_arg} '
            'VALIDATION_STRINGENCY=LENIENT '
            'REF_FLAT={input.refflat} '
            'INPUT={input.bam} '
            'OUTPUT={output.metrics} '
            '&> {log}'
        )
