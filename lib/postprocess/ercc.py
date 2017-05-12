import gzip
from snakemake.shell import shell
from Bio import SeqIO

def fasta_postprocess(origfn, newfn):
    shell(
          "gzip -c {origfn} > {newfn} "
          "&& rm {origfn}")

def gtf_postprocess(origfn, newfn):
    if isinstance(origfn, str):
        fname = origfn
    else:
        fname = origfn[0]

    with gzip.open(newfn, 'wt') as out:
        for ercc in SeqIO.parse(fname, 'fasta'):
            row = '{name}\tERCC\texon\t1\t{length}\t.\t+\t.\tgene_id "{name}"; transcript_id "{name}";\n'.format(name=ercc.id, length=len(ercc))
            out.write(row)

    shell("rm {origfn}")
