import sys
import os
import numpy as np
import pandas as pd
from snakemake.shell import shell

sys.path.insert(0, os.path.dirname(__file__) + "/..")
from lib import chipseq

# Based on the filename, identify the algorithm;
# Based on the contents, identify the format.
algorithm = os.path.basename(os.path.dirname(snakemake.input.bed))
kind = chipseq.detect_peak_format(snakemake.input.bed)

# bedToBigBed doesn't handle zero-size files
if os.stat(snakemake.input.bed).st_size == 0:
    shell("touch {output}")

# Note that autoSql filenames are relative to the workdir of the snakefile
# calling this script.
elif kind == 'narrowPeak':
    _as = '../../include/autosql/bigNarrowPeak.as'
    _type = 'bed6+4'
    names=[
        'chrom', 'chromStart', 'chromEnd', 'name', 'score',
        'strand', 'signalValue', 'pValue', 'qValue', 'peak']
elif kind == 'broadPeak':
    _as = '../../include/autosql/bigBroadPeak.as'
    _type = 'bed6+3'
    names=[
        'chrom', 'chromStart', 'chromEnd', 'name', 'score',
        'strand', 'signalValue', 'pValue', 'qValue']
elif kind == 'epic2Input':
    _as = f'../../include/autosql/{kind}Peak.as'
    _type = 'bed6+4'
    names=[
        'chrom', 'chromStart', 'chromEnd', 'pValue', 'score',
        'strand', 'ChIPCount', 'InputCount', 'FDR', 'log2FoldChange']
elif kind == 'epic2NoInput':
    _as = f'../../include/autosql/{kind}Peak.as'
    _type = 'bed6'
    names=[
        'chrom', 'chromStart', 'chromEnd', 'ChIPCount', 'score',
        'strand']
else:
    raise ValueError("Unhandled format for {0}".format(input.bed))

df = pd.read_table(snakemake.input.bed, index_col=False, names=names)
df['score'] = df['score'] - df['score'].min()
df['score'] = (df['score'] / df['score'].max()) * 1000
df['score'] = df['score'].replace([np.inf, -np.inf], np.nan).fillna(0)
df['score'] = df['score'].astype(int)
df.to_csv(snakemake.output[0] + '.tmp', sep='\t', index=False, header=False)

shell('bedToBigBed -as={_as} -type={_type} {snakemake.output}.tmp {snakemake.input.chromsizes} {snakemake.output} &> {snakemake.log}')
shell('rm {snakemake.output}.tmp')
