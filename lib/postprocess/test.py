import gzip

def test_postprocess(infiles, outfile):
    """
    converts the first chromosome name to "ABC" for testing
    """
    from lib import common
    f = common.openfile(infiles[0], 'rt')
    with gzip.open(outfile, 'wt') as fout:
        toks = f.readline().split('\t')
        toks[0] = 'ABC'
        fout.write('\t'.join(toks))
        fout.write(f.read())

