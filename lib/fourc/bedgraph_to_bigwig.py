import pybedtools
import pandas


def fixed_bedgraph(fn):
    """
    Ensure output is tab-delimited, not space-delimited.
    """
    tmp = pybedtools.BedTool._tmp()
    with open(tmp, 'w') as fout:
        for line in open(fn):
            fout.write(line.replace(' ', '\t'))
    return pybedtools.BedTool(tmp)


def weighted_mean_bedgraph(bed, bg, outfn=None, is_sorted=True):
    """
    Given a BED file and a bedGraph file with possibly overlapping regions,
    returns a new bedGraph file of the same length as the BED file which
    contains the weighted average of scores for all bedGraph regions within
    each BED region.

    If `bed` is a set of genomic windows, this is a good method for building
    bigWig-ready bedGraph files.

    Parameters
    ----------
    bed : str
        Path to BED file which will be used as "windows" into the bedGraph

    bg: str
        Path to bedGraph file.

    outfn : str or None

        TSV with fields [ chrom, start, stop, avgscore, key ] where `key` is
        "chrom:start-stop".

    is_sorted : bool
        If True, assume that the inputs are already sorted.

    Returns
    -------

    pandas.DataFrame
    """
    bg = pybedtools.BedTool(bg)
    if not is_sorted:
        bg = bg.sort()
        bed = pybedtools.BedTool(bed).sort()

    df = bg.intersect(bed, wo=True).to_dataframe(
        names=[
            'chrom',
            'start',
            'end',
            'score',
            'chrom_bed',
            'start_bed',
            'end_bed',
            'overlap',
        ])

    # group by BED region, and for each group calculate the weighted mean of
    # all possibly-overlapping bedGraph regions
    grouped = df.groupby(['chrom_bed', 'start_bed', 'end_bed'])

    if outfn is not None:
        fout = open(outfn, 'w')

    new_df = []
    for label, group in grouped:
        chrom, start, end = label
        total = float(group.overlap.sum())
        scaled_overlap = group.overlap / total
        averaged_score = (group.score * scaled_overlap).sum()
        new_df.append(
            {
                'chrom': chrom,
                'start': start,
                'end': end,
                'averaged_score': averaged_score,
            }
        )
        if outfn is not None:
            fout.write(
                '{chrom}\t{start}\t{end}\t{averaged_score}\n'.format(**locals()))
    if outfn is not None:
        fout.close()
    return pandas.DataFrame(new_df)[['chrom', 'start', 'end', 'averaged_score']]


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--windowsize', type=int)
    ap.add_argument('--bedgraph')
    ap.add_argument('--output')
    ap.add_argument('--chromsizes')
    ap.add_argument('--chromosome')
    args = ap.parse_args()

    for line in open(args.chromsizes):
        chrom, size = line.strip().split()
        if chrom == args.chromosome:
            chrom = pybedtools.BedTool(
                '{args.chromosome}\t0\t{size}\n'.format(**locals()),
                from_string=True
            )
            break

    windows = pybedtools.BedTool().window_maker(b=chrom, w=args.windowsize)
    bedgraph = fixed_bedgraph(args.bedgraph)
    weighted_mean_bedgraph(windows, bedgraph, args.output)

