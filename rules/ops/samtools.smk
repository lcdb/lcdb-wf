
rule flagstat:
    input:
        bam=patterns['markduplicates']['bam'],
        bai=patterns['markduplicates']['bam'] + '.bai'
    output:
        patterns['samtools']['flagstat']
    log:
        patterns['samtools']['flagstat'] + '.log'
    shell:
        'samtools flagstat {input.bam} > {output}'


rule samtools_stats:
    input:
        bam=patterns['markduplicates']['bam'],
        bai=patterns['markduplicates']['bam'] + '.bai'
    output:
        patterns['samtools']['stats']
    log:
        patterns['samtools']['stats'] + '.log'
    shell:
        'samtools stats {input.bam} > {output}'

rule idxstats:
    """
    Run samtools idxstats on sample bams
    """
    input:
        bam=patterns['markduplicates']['bam'],
        bai=patterns['markduplicates']['bam'] + '.bai'
    output:
        txt=patterns['samtools']['idxstats']
    log: 
        patterns['samtools']['idxstats'] + '.log'
    resources:
        mem_mb=gb(16),
        runtime=autobump(hours=2)
    run:
        shell(
            'samtools idxstats {input.bam} 2> {log} 1> {output.txt}'
        )
