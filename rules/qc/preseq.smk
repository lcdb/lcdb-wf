rule preseq:
    """
    Compute a library complexity curve with preseq
    """
    input:
        bam=patterns['bam']
    output:
        patterns['preseq']
    threads: 1
    resources:
        mem_mb=gb(1),
        runtime=autobump(hours=2)
    shell:
        'preseq '
        'c_curve '
        '-B {input} '
        '-o {output} '

