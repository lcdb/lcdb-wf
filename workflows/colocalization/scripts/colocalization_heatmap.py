import matplotlib
matplotlib.use('agg')
import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt

#outdir = snakemake.config['output']
#domain = snakemake.wildcards['domain']
#value = snakemake.wildcards['value']
#algorithm = snakemake.wildcards['algorithm']
#output = snakemake.output[0]

import argparse

ap = argparse.ArgumentParser()
ap.add_argument('--domain')
ap.add_argument('--algorithm')
ap.add_argument('--value')
ap.add_argument('--outdir')
ap.add_argument('--output')
args = ap.parse_args()

domain = args.domain
algorithm = args.algorithm
value = args.value
outdir = args.outdir
output = args.output

import sys
print('\n'.join(sys.path))

def dataframe_for_domain(domain, algorithm):
    """
    Read all files within a directory and build the dataframe.

    Empty files are listed as NaNs in the dataframe.
    """
    df = []
    files = glob.glob(os.path.join(outdir, algorithm, domain, '*', '*.txt'))
    for filename in files:
        query, reference = os.path.basename(filename).replace('.txt', '').split('_vs_')
        try:
            _df = pd.read_csv(filename, comment='#', sep='\t')
        except pd.errors.EmptyDataError:
            _df = pd.DataFrame([dict(value=np.nan)])

        _df['query'] = query
        _df['reference'] = reference
        df.append(
            _df.ix[0].to_dict()
        )
    return pd.DataFrame(df)


# Cluster methods
METRIC = 'correlation'
METHOD = 'average'


def dataframe_for_value(domain, algorithm, value):

    df = dataframe_for_domain(domain, algorithm)

    vmin, vmax = None, None

    # For IntervalStats, Use the "fraction of intervals with p<0.01" as the
    # value.
    #
    # These are all positive values. NaNs are set to 0, and the diagonal is
    # set to 1.0 (i.e., 100% of intervals are significant with respect to
    # each other)
    if algorithm == 'IntervalStats':
        piv = df.pivot(index='query', columns='reference', values=value)
        fill_piv = piv.fillna(0)
        vmax = fill_piv.max().max()
        np.fill_diagonal(fill_piv.values, 1)
        units = 'fraction pvals < 0.%s' % (value.split('_')[-1])
        title = 'IntervalStats'

    # For GAT log2foldchange, set anything with qval > 0.05 to
    # logfoldchange = 0. Diagonal is filled with 0 (log2foldchange of 1).
    # NaNs are also set to 0.
    elif algorithm == 'GAT' and value == 'l2fold':
        piv = df.pivot(index='query', columns='reference', values='l2fold')

        # used for checking
        mask = df.pivot(index='query', columns='reference',
                        values='qvalue')
        title = 'GAT foldchange'
        piv[mask > 0.05] = 0
        piv = piv.fillna(0)
        fill_piv = piv
        np.fill_diagonal(fill_piv.values, 0)
        units = 'log2fold'

    # For GAT fractions, we set the upper and lower triangles of the matrix
    # to the "track" and "annotation" overlaps in GAT terminology. We also
    # get a significance value here (qval) so we set the fraction overlap
    # to zero for anything with qval > 0.05.
    elif algorithm == 'GAT' and value == 'fractions':
        segment_frac = df.pivot(index='query', columns='reference',
                                values='percent_overlap_size_track')
        annotation_frac = df.pivot(index='query', columns='reference',
                                   values='percent_overlap_size_annotation')
        mask = df.pivot(index='query', columns='reference', values='qvalue')
        piv = segment_frac
        lower_tri_mask = np.ones(piv.shape, dtype='bool')
        lower_tri_mask[np.tril_indices(len(piv))] = False
        piv[lower_tri_mask] = annotation_frac[lower_tri_mask]
        piv[mask > 0.05] = 0
        fill = 0
        fill_piv = piv
        units = 'percentage overlap'
        title = 'GAT percentage nucleotide overlap'

    # For fisher, we want to plot the -log10(two-tail pval).
    #
    # So we keep track of the ratio, flip pvals where ratio <1, and replace
    # inf and -inf with the otherwise max and min values respectively. NaNs
    # are given a -log10(pval) = 0 (so a pval of 1.0).
    elif algorithm == 'fisher' and value == 'pval':
        piv = df.pivot(index='query', columns='reference',
                       values='two-tail')
        mask_left = df.pivot(index='query', columns='reference',
                             values='left')
        mask_right = df.pivot(index='query', columns='reference',
                              values='right')
        mask_ratio = df.pivot(index='query', columns='reference',
                              values='ratio')
        flip = mask_ratio < 1
        piv = -np.log10(piv)
        piv[flip] *= -1
        mx = piv.replace([np.inf], 0).max().max()
        mn = piv.replace([-np.inf], 0).min().min()
        piv = piv.replace([np.inf], mx)
        piv = piv.replace([-np.inf], mn)
        fill_piv = piv.fillna(0)
        units = '-log10(pval)'
        title = 'Fisher'

    ####################################################
    # TODO: also plot fisher ratio
    ####################################################

    # For jaccard, we plot the value directly. While the value can range
    # [0, 1], in practice we rarely find such good overlap.
    elif algorithm == 'jaccard' and value == 'jaccard':
        piv = df.pivot(index='query', columns='reference', values='jaccard')
        fill_piv = piv
        units = 'Jaccard statistic'
        vmin, vmax = (0, .3)
        title = 'Jaccard'

    return dict(
      fill_piv=fill_piv,
      vmin=vmin,
      vmax=vmax,
      units=units,
      title=title
    )


def plot_heatmap(fill_piv, vmin, vmax, title, units, metric='euclidean',
                 method='average', idx=None, clustermap_kwargs=dict()):
    """
    Plot a clustered heatmap of the provided values. Rows are clustered
    identically as columns so that the diagonal represents the self-self
    comparisons.

    Parameters
    ----------

    fill_piv : pandas.DataFrame
        A prepared dataframe where rownames == colnames and where -inf, inf,
        and NaN have been filled in with finite values.

    vmin, vmax : float
        Colormap limits. NOT CURRENTLY USED.

    title : str
        Title for plot

    units : str
        Units to use in colorbar

    metric : str
        Clustering metric. See `scipy.distance` for available options.

    method : clustering method
        Hierarchical clustering linkage method. See `scipy.hierarchy` for
        available options.

    idx : None or index
        If not None, then this index is used to subset `fill_piv`.

    clustermap_kwargs : dict
        Additional arguments passed to seaborn.clustermap.
    """


    fill_piv = fill_piv.astype(float)
    # subset if requested
    if idx is not None:
        fill_piv = fill_piv.ix[idx, idx]

    # Distance matrix, setting NaN to zero if necessary
    dist = distance.pdist(fill_piv.values, metric=metric)
    dist[np.isnan(dist)] = 0
    dist[dist < 0] = 0

    # ward actually uses values directly rather than using the distance matrix.
    if method == 'ward':
        vals = fill_piv.values
    else:
        vals = dist

    # Here we compute the row linkage and provide that to sns.clustermap as
    # both row and column linkages so that the same clustering is used. This
    # gets us the self-self colocalization on the diagonal.
    row_linkage = hierarchy.linkage(vals, method=method)

    # catch and fix errors in dendrogram before sending to clustermap
    mx = row_linkage[np.isfinite(row_linkage)].max()
    mn = row_linkage[np.isfinite(row_linkage)].min()
    # row_linkage[np.isinf(row_linkage)] = mx
    # scipy.clip(row_linkage, 0, mx, row_linkage)
    ind = hierarchy.dendrogram(row_linkage, no_plot=True)['leaves']


    a = sns.clustermap(fill_piv, row_linkage=row_linkage,
                       col_linkage=row_linkage, **clustermap_kwargs)

    # Fix labels
    for txt in a.ax_heatmap.get_xticklabels():
        txt.set_rotation(90)
    for txt in a.ax_heatmap.get_yticklabels():
        txt.set_rotation(0)

    # Use the provided units to label the colorbar
    a.cax.set_ylabel(units)

    # Add figure-level title and tweak margins.
    fig = plt.gcf()
    fig.suptitle(title, weight='bold', size=20)
    fig.subplots_adjust(right=0.8, bottom=0.2)
    return a


v = dataframe_for_value(domain, algorithm, value)

if (v['fill_piv'] < 0).values.any() & (v['fill_piv'] > 0).values.any():
    center = 0
    cmap = 'RdBu_r'
else:
    center = None
    cmap = sns.cubehelix_palette(as_cmap=True)


fig = plot_heatmap(
  fill_piv=v['fill_piv'],
  vmin=v['vmin'],
  vmax=v['vmax'],
  title=v['title'],
  units=v['units'],
  metric='euclidean',
  method='average',
  idx=None,
  clustermap_kwargs=dict(center=center, cmap=cmap)
)

fig.savefig(output)
