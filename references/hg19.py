from snakemake.shell import shell
def plus_lncrna_fasta_postprocess(tmpfiles, outfile):
    shell('gunzip -c {tmpfiles} > {outfile}')
