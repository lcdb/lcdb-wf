import gzip
from snakemake.shell import shell
from Bio import SeqIO

def fasta_postprocess(origfn, newfn):
    shell(
          "gzip -c {origfn} > {newfn} "
          "&& rm {origfn}")

def gtf_postprocess(origfn, newfn):
    """
    There is no GTF file for ERCC spike-ins, only a fasta. So this expects an
    input fasta file, and creates a corresponding GTF file.
    """
    if isinstance(origfn, str):
        fname = origfn
    else:
        fname = origfn[0]

    with gzip.open(newfn, 'wt') as out:
        for ercc in SeqIO.parse(fname, 'fasta'):
            row = '{name}\tERCC\texon\t1\t{length}\t.\t+\t.\tgene_id "{name}"; transcript_id "{name}";\n'.format(name=ercc.id, length=len(ercc))
            out.write(row)

    shell("rm {origfn}")


def add_gtf_to_genome(infiles, outfile, __preprocess__=None):
    processed_ercc = outfile + '.processed_ercc.tmp.gz'
    processed_non_ercc = outfile + '.processed_non_ercc.tmp.gz'

    gtf_postprocess(infiles[-1:], processed_ercc)

    if not __preprocess__:
        non_ercc = ' '.join(infiles[:-1])
        shell('cat {non_ercc} > {processed_non_ercc}')
    else:
        preprocessed_infiles = __preprocess__(infiles[:-1], processed_non_ercc)
    shell("cat {processed_non_ercc} {processed_ercc} > {outfile} "
          "&& rm {processed_non_ercc} {processed_ercc}")


def add_fasta_to_genome(infiles, outfile, reference_gzipped=True, __preprocess__=None):
    """
    Assumes that the ERCC fasta file is the *last* one provided.

    Parameters
    ----------

    infiles : str or list
        Input files to process. The last file is assumed to be the ERCC file.

    outfile : str
        Gzipped output file to process

    preprocess : callable
        Function which will be called 
    """
    non_ercc = ' '.join(infiles[:-1])
    ercc = infiles[-1]
    if reference_gzipped:
        cmd = (
            f"cat {non_ercc} > {outfile} "
            f"&& gzip -c {ercc} >> {outfile} "
            f"&& rm {non_ercc} {ercc}"
        )
        shell(cmd)
    else:
        shell("cat {infiles} | gzip -c > {outfile} && rm {infiles}")
