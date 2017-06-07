import pybedtools
from matplotlib import pyplot as plt
import matplotlib
from pybedtools import featurefuncs as ff
import numpy as np
import pandas
import bedgraph_to_bigwig as b2w
from collections import defaultdict

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('--bed')
ap.add_argument('--output')
ap.add_argument('--control', nargs='+')
ap.add_argument('--treatment', nargs='+')
args = ap.parse_args()


bed = args.bed
output = args.output
bgs = {
    'control': sorted(args.control),
    'treatment': sorted(args.treatment),
}

lookup = {}
inc = defaultdict(int)
for k, vs in bgs.items():
    for v in vs:
        num = inc[k]
        inc[k] += 1
        print(v)
        df = b2w.weighted_mean_bedgraph(bed=bed, bg=b2w.fixed_bedgraph(v))
        if 'chrom' not in lookup:
            lookup['chrom'] = df['chrom']
            lookup['start'] = df['start']
            lookup['end'] = df['end']
        lookup[k + str(num)] = df['averaged_score']

agg = pandas.DataFrame(lookup)
control_cols = [i for i in lookup.keys() if 'control' in i]
treatment_cols = [i for i in lookup.keys() if 'treatment' in i]

agg['fc'] = agg[treatment_cols].sum(axis=1) / agg[control_cols].sum(axis=1)
agg['lfc'] = np.log2(agg['fc'])
agg['name'] = '.'
tmp = pybedtools.BedTool._tmp()
agg[['chrom', 'start', 'end', 'name', 'lfc']].to_csv(
    tmp, sep='\t', index=False, header=False)

bt = pybedtools.BedTool(tmp)
norm = bt.colormap_normalize()
norm.vmax = 3
norm.vmin = -3

def fix(x):
  lfc = x.score
  return pybedtools.create_interval_from_list(list(map(str, [
    x.chrom,
    x.start,
    x.stop,
    '.',
    '0',
    '.',
    x.start,
    x.stop,
    x[8],
    lfc])))

bt = (
  bt
  .each(ff.add_color, cmap=matplotlib.cm.RdBu_r, norm=norm)
  .each(fix)
  .saveas(output)
)

