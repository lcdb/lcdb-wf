import gzip
import glob
import tarfile
from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    """
    The fasta from UCSC comes as a tarball of fastas. So we extract them all to
    a temp directory and then cat them all together into the final fa.gz file.
    """
    assert (
        (isinstance(origfn, list)) and (len(origfn) == 1)
    ), 'unexpected input: %s' % origfn
    origfn = origfn[0]
    t = tarfile.open(origfn)
    shell('mkdir -p {origfn}.tmp')
    t.extractall(origfn + '.tmp')
    with gzip.open(newfn, 'wt') as fout:
        for fa in sorted(glob.glob(origfn + '.tmp/*.fa')):
            print(fa)
            fout.write(open(fa).read())
    shell('rm -r {origfn}.tmp')
