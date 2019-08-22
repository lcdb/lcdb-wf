import gzip

def test_postprocess(infiles, outfile):
    """
    converts the first feature's featuretype to "ABC" for testing
    postprocessing
    """
    from lib import common
    f = common.openfile(infiles[0], 'rt')
    with gzip.open(outfile, 'wt') as fout:
        toks = f.readline().split('\t')
        toks[2] = 'ABC'
        fout.write('\t'.join(toks))
        fout.write(f.read())

