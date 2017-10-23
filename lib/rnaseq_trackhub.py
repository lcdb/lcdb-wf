#!/usr/bin/env python

"""
Build and upload a track hub of RNA-seq signals, with samples, sort order,
selection matrix, colors, and server info configured in a hub config file.

This module assumes a particular filename pattern and whether or not bigwigs
are stranded.  Such assumptions are indicated in the comments below.
"""

import re
import os
import pandas
import yaml
import matplotlib
from trackhub.helpers import sanitize
from trackhub import CompositeTrack, ViewTrack, SubGroupDefinition, Track, default_hub
from trackhub.upload import upload_hub

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('config', help='Main config.yaml file')
args = ap.parse_args()

# Access configured options. See comments in example hub_config.yaml for
# details
config = yaml.load(open(args.config))
hub_config_fn = os.path.join(os.path.dirname(args.config), config['hub_config'])
hub_config = yaml.load(open(hub_config_fn))


hub, genomes_file, genome, trackdb = default_hub(
    hub_name=hub_config['hub']['name'],
    short_label=hub_config['hub']['short_label'],
    long_label=hub_config['hub']['long_label'],
    email=hub_config['hub']['email'],
    genome=hub_config['hub']['genome']
)

#hub.url = hub_config['hub']['url']
#hub.remote_fn = hub_config['hub']['remote_fn']

# Set up subgroups based on the configured columns
df = pandas.read_table(config['sampletable'], comment='#')
cols = hub_config['subgroups']['columns']
subgroups = []
for col in cols:
    unique = list(df[col].unique())
    s = SubGroupDefinition(
        name=sanitize(col, strict=True),
        label=col,
        mapping={
            sanitize(i, strict=True): sanitize(i, strict=False)
            for i in unique}
    )
    subgroups.append(s)

# also add direction as a subgroup
# ASSUMPTION: stranded bigwigs were created
subgroups.append(
    SubGroupDefinition(
        name='strand',
        label='strand',
        mapping={'pos': 'pos', 'neg': 'neg'}))


def dimensions_from_subgroups(s):
    """
    Given a sorted list of subgroups, return a string appropriate to provide as
    a composite track's `dimensions` arg
    """
    letters = 'XYABCDEFGHIJKLMNOPQRSTUVWZ'
    return ' '.join(['dim{0}={1}'.format(dim, sg.name) for dim, sg in zip(letters, s)])


def filter_composite_from_subgroups(s):
    """
    Given a sorted list of subgroups, return a string appropriate to provide as
    the a composite track's `filterComposite` argumen argumen
    """
    dims = []
    for letter, sg in zip('ABCDEFGHIJKLMNOPQRSTUVWZ', s[2:]):
        dims.append('dim{0}'.format(letter))
    if dims:
        return ' '.join(dims)

# Identify the sort order based on the config, and create a string appropriate
# for use as the `sortOrder` argument of a composite track.
to_sort = hub_config['subgroups'].get('sort_order', [])
to_sort += [sg.name for sg in subgroups if sg.name not in to_sort]
sort_order = ' '.join([i + '=+' for i in to_sort])

# Identify samples based on config and sampletable
sample_dir = config['sample_dir']
samples = df[df.columns[0]]

composite = CompositeTrack(
    name=hub_config['hub']['name'] + 'composite',
    short_label='rnaseq composite',
    long_label='rnaseq composite',
    dimensions=dimensions_from_subgroups(subgroups),
    filterComposite=filter_composite_from_subgroups(subgroups),
    sortOrder=sort_order,
    tracktype='bigWig')

# ASSUMPTION: stranded bigwigs
pos_signal_view = ViewTrack(
    name='possignalviewtrack', view='possignal', visibility='full',
    tracktype='bigWig', short_label='plus strand', long_label='plus strand signal')
neg_signal_view = ViewTrack(
    name='negsignalviewtrack', view='negsignal', visibility='full',
    tracktype='bigWig', short_label='minus strand', long_label='minus strand signal')

supplemental_view = ViewTrack(
    name='suppviewtrack', view='supplementa', visibility='full',
    tracktype='bigBed', short_label='Supplemental', long_label='Supplemental')

colors = hub_config.get('colors', [])


def hex2rgb(h):
    """
    Given a hex color code, return a 0-255 RGB tuple as a CSV string, e.g.,
    "#ff0000" -> "255,0,0"
    """
    return ','.join(map(lambda x: str(int(x * 255)), matplotlib.colors.hex2color(h)))


def decide_color(samplename):
    """
    Look up the color dictionary in the config and return the first color that
    matches `samplename`.
    """
    for cdict in hub_config.get('colors', []):
        k = list(cdict.keys())
        assert len(k) == 1
        color = k[0]
        v = cdict[color]
        for pattern in v:
            regex = re.compile(pattern)
            if regex.search(samplename):
                return hex2rgb(color)
    return hex2rgb('#000000')


for sample in df[df.columns[0]]:
    # ASSUMPTION: stranded bigwigs
    for direction in 'pos', 'neg':

        # ASSUMPTION: bigwig filename pattern
        bigwig = os.path.join(
            sample_dir, sample,
            sample + '.cutadapt.bam.{0}.bigwig'.format(direction))

        subgroup = df[df.iloc[:, 0] == sample].to_dict('records')[0]
        subgroup = {
            sanitize(k, strict=True): sanitize(v, strict=True)
            for k, v in subgroup.items()
        }

        # ASSUMPTION: stranded bigwigs
        additional_kwargs = {}
        subgroup['strand'] = direction
        view = pos_signal_view
        if direction == 'neg':
            # additional_kwargs['negateValues'] = 'on'
            # additional_kwargs['viewLimits'] = '-25:0'
            additional_kwargs['viewLimits'] = '15:0'
            view = neg_signal_view
        else:
            # if strands were switched....
            additional_kwargs['negateValues'] = 'on'
            additional_kwargs['viewLimits'] = '-15:0'
            # additional_kwargs['viewLimits'] = '0:25'
        view.add_tracks(
            Track(
                name=sanitize(sample + os.path.basename(bigwig), strict=True),
                short_label=sample + '_' + direction,
                long_label=sample + '_' + direction,
                tracktype='bigWig',
                subgroups=subgroup,
                source=bigwig,
                color=decide_color(sample),
                altColor=decide_color(sample),
                maxHeightPixels='8:35:100',
                **additional_kwargs
            )
        )

supplemental = hub_config.get('supplemental', [])
if supplemental:
    composite.add_view(supplemental_view)
    for block in supplemental:
        supplemental_view.add_tracks(Track(**block))


# Tie everything together
composite.add_subgroups(subgroups)
trackdb.add_tracks(composite)
composite.add_view(pos_signal_view)
composite.add_view(neg_signal_view)

# Render and upload using settings from hub config file
hub.render()
kwargs = hub_config.get('upload', {})
upload_hub(hub=hub, **kwargs)
print(hub.url)
