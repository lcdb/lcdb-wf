from Bio import SeqIO
import gzip
from snakemake.shell import shell

def rrna_postprocess(tmpfiles, outfile):
    def gen():
        for tmp in tmpfiles:
            handle = gzip.open(tmp, 'rt')
            parser = SeqIO.parse(handle, 'fasta')
            for rec in parser:
                if 'Dictyostelium discoideum' not in rec.description:
                    continue
                rec.seq = rec.seq.back_transcribe()
                rec.description = rec.name
                yield rec

    with gzip.open(outfile, 'wt') as fout:
        SeqIO.write(gen(), fout, 'fasta')
