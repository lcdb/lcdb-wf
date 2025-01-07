# TODO: split into multiple featurecounts runs, since PE needs to be sorted each time.
rule featurecounts:
    """
    Count reads in annotations with featureCounts from the subread package
    """
    input:
        annotation=rules.gtf.output,
        bam=expand(patterns['markduplicates']['bam'], sample=SAMPLES),
    output:
        counts=patterns['featurecounts']
    log:
        patterns['featurecounts'] + '.log'
    threads: 8
    resources:
        mem_mb=gb(16),
        runtime=autobump(hours=2)
    params:
        strand_arg={
                'unstranded': '-s0 ',
                'fr-firststrand': '-s2 ',
                'fr-secondstrand': '-s1 ',
            }[config["stranded"]],
        extra=""
    run:
        # NOTE: By default, we use -p for paired-end
        p_arg = ''
        if is_paired:
            p_arg = '-p --countReadPairs '
        shell(
            'featureCounts '
            '{params.strand_arg} '
            '{p_arg} '
            '-T {threads} '
            '-a {input.annotation} '
            '-o {output.counts} '
            '{input.bam} '
            '&> {log}'
        )
