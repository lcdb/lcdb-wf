
rule bigwig_neg:
    """
    Create a bigwig for negative-strand reads
    """
    input:
        bam=patterns['markduplicates']['bam'],
        bai=patterns['markduplicates']['bam'] + '.bai',
    output:
        patterns['bigwig']['neg']
    threads: 8
    resources:
        mem_mb=gb(16),
        runtime=autobump(hours=2)
    log:
        patterns['bigwig']['neg'] + '.log'
    params:
        strand_arg = {
                'unstranded': '',
                'fr-firststrand': '--filterRNAstrand reverse ',
                'fr-secondstrand': '--filterRNAstrand forward ',
            }[config["stranded"]],
        extra=(
            '--minMappingQuality 20 '
            '--smoothLength 10 '
            '--normalizeUsing BPM '    # equivalent to TPM # [TEST SETTINGS]
        ),
    run:
        shell(
            'bamCoverage '
            '--bam {input.bam} '
            '-o {output} '
            '-p {threads} '
            '{params.extra} '
            '{params.strand_arg} '
            '&> {log}'
        )


rule bigwig_pos:
    """
    Create a bigwig for postive-strand reads.
    """
    input:
        bam=patterns['markduplicates']['bam'],
        bai=patterns['markduplicates']['bam'] + '.bai',
    output:
        patterns['bigwig']['pos']
    threads: 8
    resources:
        mem_mb=gb(16),
        runtime=autobump(hours=2)
    log:
        patterns['bigwig']['pos'] + '.log'
    params:
        strand_arg={
                'unstranded': '',
                'fr-firststrand': '--filterRNAstrand forward ',
                'fr-secondstrand': '--filterRNAstrand reverse ',
            }[config["stranded"]],
        extra=(
            '--minMappingQuality 20 '
            '--smoothLength 10 '
            '--normalizeUsing BPM '    # equivalent to TPM # [TEST SETTINGS]
        ),
    run:
        shell(
            'bamCoverage '
            '--bam {input.bam} '
            '-o {output} '
            '-p {threads} '
            '{BAMCOVERAGE_ARGS} '
            '{params.strand_arg} '
            '&> {log}'
        )

