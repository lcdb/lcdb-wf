
# STAR can be run in 1-pass or 2-pass modes. Since we may be running it
# more than once in almost the same way, we pull out the shell command here
# and use it below.
STAR_CMD = (
    'STAR '
    '--runThreadN {threads} '
    '--genomeDir {genomedir} '
    '--readFilesIn {input.fastq} '
    '--readFilesCommand zcat '
    '--outFileNamePrefix {prefix} '
    '{params.extra} '
)
STAR_PARAMS = (
    # NOTE: The STAR docs indicate that the following parameters are
    # standard options for ENCODE long-RNA-seq pipeline.  Comments are from
    # the STAR docs.
    '--outFilterType BySJout '               # reduces number of spurious junctions
    '--outFilterMultimapNmax 20 '            # if more than this many multimappers, consider unmapped
    '--alignSJoverhangMin 8 '                # min overhang for unannotated junctions
    '--alignSJDBoverhangMin 1 '              # min overhang for annotated junctions
    '--outFilterMismatchNmax 999 '           # max mismatches per pair
    '--outFilterMismatchNoverReadLmax 0.04 ' # max mismatches per pair relative to read length
    '--alignIntronMin 20 '                   # min intron length
    '--alignIntronMax 1000000 '              # max intron length
    '--alignMatesGapMax 1000000 '            # max distance between mates
    '--outSAMunmapped None '                 # do not report aligned reads in output
)
logfile_extensions = ['Log.progress.out', 'Log.out', 'Log.final.out', 'Log.std.out']

if config['aligner'] == 'star':

    rule star:
        """
        Align with STAR (1-pass mode)
        """
        input:
            fastq=expand(patterns["cutadapt"], n=n, allow_missing=True),
            index=rules.star_index.output,
            annotation=f"{REFERENCES}/annotation.gtf"
        output:
            bam=temporary(patterns['bam']),
            sjout=temporary(patterns['bam'].replace('.bam', '.star.SJ.out.tab')),
        log:
            patterns['bam'].replace('.bam', '.star.bam.log')
        threads: 16
        resources:
            mem_mb=gb(64),
            runtime=autobump(hours=8)
        params:
            extra=STAR_PARAMS

        run:
            genomedir = os.path.dirname(input.index[0])
            outdir = os.path.dirname(output[0])
            prefix = output.bam.replace('.bam', '.star.')
            shell(
                STAR_CMD + (
                    '--outSAMtype BAM SortedByCoordinate '
                    '--outStd BAM_SortedByCoordinate > {output.bam} '
                    '2> {log} '
                )
            )

            # move various hard-coded log files to log directory
            logfiles = expand(prefix + '{ext}', ext=logfile_extensions)
            shell('mkdir -p {outdir}/star_logs '
                  '&& mv {logfiles} {outdir}/star_logs')

if config['aligner'] == 'star-twopass':

    rule star_pass1:
        """
        First pass of alignment with STAR to get the junctions
        """
        input:
            fastq=expand(patterns["cutadapt"], n=n, allow_missing=True),
            index=rules.star_index.output,
            annotation=f"{REFERENCES}/annotation.gtf"
        output:
            sjout=temporary(patterns['bam'].replace('.bam', '.star-pass1.SJ.out.tab')),
        log:
            patterns['bam'].replace('.bam', '.star-pass1.bam.log')
        threads: 16
        resources:
            mem_mb=gb(64),
            runtime=autobump(hours=8)
        params:
            extra=STAR_PARAMS
        run:
            genomedir = os.path.dirname(input.index[0])
            outdir = os.path.dirname(output[0])
            prefix = output.sjout.replace('SJ.out.tab', '')
            shell(
                STAR_CMD +
                (
                    # In this first pass, we don't actually care about the
                    # alignment -- just the detected junctions. So we output
                    # the SAM to /dev/null.
                    '--outStd SAM > /dev/null '
                    '2> {log} '
                )
            )

            # move various hard-coded log files to log directory
            logfiles = expand(prefix + '{ext}', ext=logfile_extensions)
            shell('mkdir -p {outdir}/star-pass1_logs '
                  '&& mv {logfiles} {outdir}/star-pass1_logs')


    rule star_pass2:
        """
        Second pass of alignment with STAR using splice junctions across all
        samples to get the final BAM
        """
        input:
            fastq=expand(patterns["cutadapt"], n=n, allow_missing=True),
            index=rules.star_index.output,
            annotation=f"{REFERENCES}/annotation.gtf",
            sjout=expand(patterns['bam'].replace('.bam', '.star-pass1.SJ.out.tab'), sample=SAMPLES),
        output:
            bam=temporary(patterns['bam']),
            sjout=temporary(patterns['bam'].replace('.bam', '.star-pass2.SJ.out.tab')),
        log:
            patterns['bam'].replace('.bam', '.star-pass2.bam.log')
        threads: 16
        resources:
            mem_mb=gb(64),
            runtime=autobump(hours=8)
        params:
            extra=STAR_PARAMS
        run:
            genomedir = os.path.dirname(input.index[0])
            outdir = os.path.dirname(output[0])
            prefix = output.bam.replace('.bam', '.star-pass2.')
            shell(
                STAR_CMD + (
                    # In contrast to pass 1, we will be keeping these BAMs --
                    # so sort them
                    '--outSAMtype BAM SortedByCoordinate '

                    # Splice junction databases from all samples in the first
                    # pass.
                    '--sjdbFileChrStartEnd {input.sjout} '
                    '--outStd BAM_SortedByCoordinate > {output.bam} '
                    '2> {log} '
                )
            )

            # move various hard-coded log files to log directory
            logfiles = expand(prefix + '{ext}', ext=logfile_extensions)
            shell('mkdir -p {outdir}/star-pass2_logs '
                  '&& mv {logfiles} {outdir}/star-pass2_logs')

            shell('rm -r {prefix}_STARgenome')
