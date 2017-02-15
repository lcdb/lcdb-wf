from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    shell(
          "unzip -p {origfn} -x ERCC92.gtf | "
          "gzip -c > {newfn} "
          "&& rm {origfn}")

def gtf_postprocess(origfn, newfn):
    shell(
          "unzip -p {origfn} -x ERCC92.fa | "
          "gzip -c > {newfn} "
          "&& rm {origfn}")
