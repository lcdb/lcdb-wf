from snakemake.shell import shell
import os

def fasta_postprocess(origfn, newfn):
    tmpdir = os.path.join(os.path.dirname(newfn), 'tmp')
    assert len(origfn) == 1
    origfn = origfn[0]
    rel_orig = os.path.relpath(origfn, tmpdir)
    rel_new = os.path.relpath(newfn, tmpdir)
    shell(
        "mkdir -p {tmpdir} "
        "&& cd {tmpdir} "
        "&& tar -xf {rel_orig} --no-same-owner "
        "&& cat PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa | "
        "fold -w 80 | "
        "gzip -c > {rel_new} "
        "&& rm {rel_orig} "
        "&& cd - && rm -r {tmpdir} "
    )
