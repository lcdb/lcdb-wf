from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    shell(
          "gzip -c {origfn} > {newfn} "
          "&& rm {origfn}")
