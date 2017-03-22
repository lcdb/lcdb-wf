from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    shell(
          "tar -xf {origfn} "
          "&& cat PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa | "
          "fold -w 80 | "
          "gzip -c > {newfn} "
          "&& rm {origfn}")
