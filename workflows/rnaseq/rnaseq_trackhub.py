#!/usr/bin/env python

"""
Build and upload a track hub of RNA-seq signals, with samples, sort order,
selection matrix, colors, and server info configured in a hub config file.

This module assumes a particular filename pattern and whether or not bigwigs
are stranded.  Such assumptions are indicated in the comments below.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import re
from pprint import pprint
import pandas
import yaml
import matplotlib
from snakemake.utils import update_config
from trackhub.helpers import sanitize, hex2rgb, dimensions_from_subgroups, filter_composite_from_subgroups
from trackhub import CompositeTrack, ViewTrack, SubGroupDefinition, Track, default_hub
from trackhub.upload import upload_hub, stage_hub
import argparse

from lib.patterns_targets import RNASeqConfig

ap = argparse.ArgumentParser()
ap.add_argument('config', help='Main config.yaml file')
ap.add_argument('hub_config', help='Track hub config YAML file')
ap.add_argument('--additional-configs', nargs='+',
                help='Additional config files with which to update the main '
                'config',)
args = ap.parse_args()

# Access configured options. See comments in example hub_config.yaml for
# details
config = yaml.load(open(args.config), Loader=yaml.FullLoader)
hub_config = yaml.load(open(args.hub_config), Loader=yaml.FullLoader)

if args.additional_configs:
    for cfg in args.additional_configs:
        update_config(config, yaml.load(open(cfg), Loader=yaml.FullLoader))

c = RNASeqConfig(config, os.path.join(os.path.dirname(args.config), 'rnaseq_patterns.yaml'))

hub, genomes_file, genome, trackdb = default_hub(
    hub_name=hub_config['hub']['name'],
    short_label=hub_config['hub']['short_label'],
    long_label=hub_config['hub']['long_label'],
    email=hub_config['hub']['email'],
    genome=hub_config['hub']['genome']
)

# Set up subgroups based on the configured columns
df = pandas.read_csv(config['sampletable'], comment='#', sep='\t')
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


# Identify the sort order based on the config, and create a string appropriate
# for use as the `sortOrder` argument of a composite track.
to_sort = hub_config['subgroups'].get('sort_order', [])
to_sort += [sg.name for sg in subgroups if sg.name not in to_sort]
sort_order = ' '.join([i + '=+' for i in to_sort])

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
    tracktype='bigWig', short_label='pos strand', long_label='positive strand signal')
neg_signal_view = ViewTrack(
    name='negsignalviewtrack', view='negsignal', visibility='full',
    tracktype='bigWig', short_label= 'neg strand', long_label='negative strand signal')

supplemental_view = ViewTrack(
    name='suppviewtrack', view='supplementalview', visibility='full',
    tracktype='bigBed', short_label='Supplemental', long_label='Supplemental')

colors = hub_config.get('colors', [])


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
        bigwig = c.patterns['bigwig'][direction].format(sample=sample)

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
            additional_kwargs['negateValues'] = 'on'
            additional_kwargs['viewLimits'] = '-25:0'
            view = neg_signal_view
        else:
            additional_kwargs['viewLimits'] = '0:25'
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

# If the hub config specified a remote path, upload there -- otherwise don't do
# any uploading.
if kwargs.get('remote_dir', False):
    upload_hub(hub=hub, **kwargs)
else:
    stage_hub(hub, 'staging')
