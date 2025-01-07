rule salmon:
    """
    Quantify reads coming from transcripts with Salmon
    """
    input:
        fastq=expand(patterns["cutadapt"], n=n, allow_missing=True),
        index=REFERENCES + "/salmon/versionInfo.json"
    output:
        patterns['salmon']
    log:
        patterns['salmon'] + '.log'
    params:
        extra=(
            "--libType=A "
            "--gcBias "
            "--seqBias "
            "--validateMappings "
        )
    threads: 6
    resources:
        mem_mb=gb(32),
        runtime=autobump(hours=2)
    run:
        outdir = os.path.dirname(output[0])
        index_dir = os.path.dirname(input.index)
        if is_paired:
            fastq_arg = f'-1 {input.fastq[0]} -2 {input.fastq[1]} '
        else:
            fastq_arg = f'-r {input.fastq} '
        shell(
            'salmon quant '
            '--index {index_dir} '
            '--output {outdir} '
            '--threads {threads} '
            '{params.extra} '
            '{fastq_arg} '
            '&> {log}'
        )
