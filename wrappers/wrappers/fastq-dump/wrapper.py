from snakemake import shell
output = snakemake.output
log = snakemake.log

srr = snakemake.params.sampletable.loc[snakemake.wildcards.sample, 'Run']

if hasattr(snakemake.params, "limit"):
    limit = f'-X {snakemake.params.limit}'
else:
    limit = ""

# Two different paths depending on the layout. In both cases, we
# want to avoid creating the final output until the very end, to
# avoid incomplete downloads.
if snakemake.params.is_paired:
    # For PE we need to use --split-files, which also means using
    # the slower --gzip
    shell(
        'fastq-dump '
        '{srr} '
        '--gzip '
        '--split-files '
        '{limit} '
        '&> {log}'
    )

    # The filenames are predictable, so we can move them as needed.
    shell('mv {srr}_1.fastq.gz {output[0]}')
    shell('mv {srr}_2.fastq.gz {output[1]}')

else:
    # For SE, we can use the faster stdout | gzip, and move it
    # directly when done.
    shell(
        'fastq-dump '
        '{srr} '
        '-Z '
        '{limit} '
        '2> {log} | gzip -c > {output[0]}.tmp '
        '&& mv {output[0]}.tmp {output[0]} '
    )
