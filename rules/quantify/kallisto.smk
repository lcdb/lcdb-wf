rule kallisto:
    """
    Quantify reads coming from transcripts with Kallisto
    """
    input:
        fastq=expand(patterns["cutadapt"], n=n, allow_missing=True),
        index=REFERENCES + "/kallisto/transcripts.idx",
    output:
        patterns['kallisto']
    params:
        strand_arg={
                'unstranded': '',
                'fr-firststrand': '--rf-stranded',
                'fr-secondstrand': '--fr-stranded',
            }[config["stranded"]],
        extra=(
            "--bootstrap-samples 100" if is_paired else
            "--single --fragment-length 300 --sd 20 --bootstrap-samples 100"
        ),
    log:
        patterns['kallisto'] + '.log'
    threads:
        8
    resources:
        mem_mb=gb(32),
        runtime=autobump(hours=2),
    run:
        outdir = os.path.dirname(output[0])
        shell(
            'kallisto quant '
            '--index {input.index} '
            '--output-dir {outdir} '
            '--threads {threads} '
            '--bootstrap-samples 100 '
            '--threads {threads} '
            '{params.strand_arg} '
            '{params.extra} '
            '{input.fastq} '
            '&> {log}'
        )
