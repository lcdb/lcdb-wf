import pybedtools
import gzip
from snakemake.shell import shell
import os


def strip_ensembl_version(infiles, outfile):
    def transform(f):
        f.attrs['gene_id'] = f.attrs['gene_id'].split('.')[0]
        return f
    with gzip.open(outfile, 'wt') as fout:
        for infile in infiles:
            for feature in pybedtools.BedTool(infile):
                fout.write(str(transform(feature)))
