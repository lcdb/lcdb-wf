
rule rseqc_infer_experiment:
    """
    Infer strandedness of experiment
    """
    input:
        bam=patterns['markduplicates']['bam'],
        bed12=rules.conversion_bed12.output,
    output:
        txt=patterns['rseqc']['infer_experiment']
    log:
        patterns['rseqc']['infer_experiment'] + '.log'
    resources:
        mem_mb=gb(2),
        runtime=autobump(hours=2)
    shell:
        'infer_experiment.py -r {input.bed12} -i {input.bam} > {output} &> {log}'


rule rseqc_read_distribution:
    """
    read distribution plots
    """
    input:
        bam=patterns['markduplicates']['bam'],
        bed12=rules.conversion_bed12.output,
    output:
        txt=patterns['rseqc']['read_distribution']
    log:
        patterns['rseqc']['read_distribution'] + '.log'
    resources:
        mem_mb=gb(2),
        runtime=autobump(hours=2)
    shell:
        'read_distribution.py -i {input.bam} -r {input.bed12} > {output} &> {log}'
