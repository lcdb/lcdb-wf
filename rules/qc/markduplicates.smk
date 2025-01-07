rule markduplicates:
    """
    Mark or remove PCR duplicates with Picard MarkDuplicates
    """
    input:
        bam=patterns['bam']
    output:
        bam=patterns['markduplicates']['bam'],
        metrics=patterns['markduplicates']['metrics'],
    log:
        patterns['markduplicates']['bam'] + '.log'
    params:
        java_args='-Xmx20g'
        # java_args='-Xmx2g'  # [TEST SETTINGS -1]
    threads: 1
    resources:
        mem_mb=gb(32),
        runtime=autobump(hours=2),
        disk_mb=autobump(gb=100),
    shell:
        'picard '
        '{params.java_args} '
        'MarkDuplicates '
        'INPUT={input.bam} '
        'OUTPUT={output.bam} '
        'METRICS_FILE={output.metrics} '
        'VALIDATION_STRINGENCY=LENIENT '
        '&> {log}'
