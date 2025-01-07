rule rRNA:
    input:
        fastq=expand(patterns["cutadapt"], n=1, allow_missing=True),  # currently only R1
        index=multiext(
            f"{REFERENCES}/bowtie2/rrna",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
            ".fa",
        ),
    output:
        bam=temporary(patterns['rrna']['bam'])
    log:
        patterns['rrna']['bam'] + '.log'
    threads: 6
    resources:
        mem_mb=gb(2),
        runtime=autobump(hours=2)
    params:
        extra=(
            '-k 1 '       # NOTE: we only care if >=1 mapped
            '--no-unal '  # NOTE: suppress unaligned reads
        )
    run:
        prefix = os.path.commonprefix(input.index).rstrip(".")
        sam = output.bam.replace('.bam', '.sam')

        shell(
            "bowtie2 "
            "-x {prefix} "
            "-U {input.fastq} "
            "--threads {threads} "
            "{params.extra} "
            "-S {sam} "
            "> {log} 2>&1"
        )

        shell(
            "samtools view -Sb {sam} "
            "| samtools sort - -o {output.bam} -O BAM "
            "&& rm {sam}"
        )

rule rrna_libsizes_table:
    """
    Aggregate rRNA counts into a table
    """
    input:
        rrna=expand(patterns['rrna']['libsize'], sample=SAMPLES),
        fastq=expand(patterns['libsizes']['cutadapt'], sample=SAMPLES),
    output:
        json=patterns['rrna_percentages_yaml'],
        tsv=patterns['rrna_percentages_table']
    threads: 1
    resources:
        mem_mb=gb(2),
        runtime=autobump(hours=2)
    run:
        def rrna_sample(f):
            return utils.extract_wildcards(patterns['rrna']['libsize'], f)['sample']

        def sample(f):
            return utils.extract_wildcards(patterns['libsizes']['cutadapt'], f)['sample']

        def million(f):
            return float(open(f).read()) / 1e6

        rrna = sorted(input.rrna, key=rrna_sample)
        fastq = sorted(input.fastq, key=sample)
        samples = list(map(rrna_sample, rrna))
        rrna_m = list(map(million, rrna))
        fastq_m = list(map(million, fastq))

        df = pd.DataFrame(dict(
            sample=samples,
            million_reads_rRNA=rrna_m,
            million_reads_fastq=fastq_m,
        ))
        df = df.set_index('sample')
        df['rRNA_percentage'] = df.million_reads_rRNA / df.million_reads_fastq * 100

        df[['million_reads_fastq', 'million_reads_rRNA', 'rRNA_percentage']].to_csv(output.tsv, sep='\t')
        y = {
            'id': 'rrna_percentages_table',
            'section_name': 'rRNA content',
            'description': 'Amount of reads mapping to rRNA sequence',
            'plot_type': 'table',
            'pconfig': {
                'id': 'rrna_percentages_table_table',
                'title': 'rRNA content table',
                'min': 0
            },
            'data': yaml.load(df.transpose().to_json(), Loader=yaml.FullLoader),
        }
        with open(output.json, 'w') as fout:
            yaml.dump(y, fout, default_flow_style=False)
