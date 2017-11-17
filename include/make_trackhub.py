import os
import argparse
from collections import defaultdict
import re
import yaml
from trackhub import default_hub
from trackhub import upload
from trackhub import Track, CompositeTrack, ViewTrack, AggregateTrack, SuperTrack
from trackhub.track import SubGroupDefinition
from trackhub.helpers import sanitize


def hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

ap = argparse.ArgumentParser()
ap.add_argument('--config')
args = ap.parse_args()


# figure out from the config what files we will have
config = yaml.load(open(args.config))

colors = config['trackhub']['colors']
color_regexp = []
for i in colors:
    for k, v in i.items():
        color_regexp.append(
            (re.compile(k), hex_to_rgb(v))
        )

def color_for_track(fn):
    for regexp, color in color_regexp:
        if regexp.search(fn):
            return ','.join(list(map(str, color)))
    return '0,0,0'

hub, genomes_file, genome, trackdb = default_hub(**config['trackhub']['hub'])

hub.remote_fn = config['trackhub']['upload']['remote_fn']

sample_composite = CompositeTrack(
    name='samplespecific', short_label='Samples',
    long_label='Individual samples', tracktype='bigWig')

sample_signal_view = ViewTrack(
    name='samplesignalview', view='samplesignal', visibility='full', tracktype='bigWig',
    viewLimits='0:150000', maxHeightPixels='15:50:127',
    short_label='Sample Signal', long_label='Sample-level Signal')

sample_bed_view = ViewTrack(
    name='samplebedview', view='samplebed', visibility='dense', tracktype='bigBed',
    short_label='Sample Region', long_label='Sample-level regions')


bait_composite = CompositeTrack(
    name='baitspecific', short_label='Baits', long_label='Bait level',
    tracktype='bigBed')

bait_bed_view = ViewTrack(name='bedview', view='bed', visibility='squish',
                          tracktype='bigBed 3', short_label='Regions',
                          long_label='Regions')

supertrack = SuperTrack(
    name='supertrack',
    short_label='Aggregated sample supertrack',
    long_label='Aggregated sample supertrack',
    viewLimits='0:150000',

)


sample_composite.add_view(sample_signal_view)
sample_composite.add_view(sample_bed_view)

bait_composite.add_view(bait_bed_view)

trackdb.add_tracks(sample_composite)
trackdb.add_tracks(bait_composite)
trackdb.add_tracks(supertrack)

# We will build the possible subgroups up from the track metadata, and build
# SubGroupDefinitions later.
sample_subgroup_values = defaultdict(set)
bait_subgroup_values = defaultdict(set)

aggregates = {}

for comparison, vals, in config['4c']['comparisons'].items():
    for kind in ['cis', 'nearbait']:

        bait = vals['bait']
        k = config['4c']['baits'][bait][kind + '_k']

        for treatment in ['control', 'treatment', 'all']:


            bait_subgroups = dict(
                bait=sanitize(bait),
                treatment=sanitize(treatment),
                kind=sanitize(kind),
            )
            for key, val in bait_subgroups.items():
                bait_subgroup_values[sanitize(key)].update(
                    [sanitize(val, False)])

            if treatment == 'all':
                local_fn = (
                    '4cker-output/{comparison}/{kind}_k{k}/'
                    '{bait}_{kind}_k{k}_adaptive_windows.bigbed'.format(**locals())
                )

                bait_bed_view.add_tracks(
                    Track(
                        name=sanitize(os.path.basename(local_fn)),
                        tracktype='bigBed 3',
                        local_fn=local_fn,
                        short_label=os.path.basename(local_fn),
                        long_label=os.path.basename(local_fn),
                        subgroups=bait_subgroups,
                        color=color_for_track(local_fn),
                    )
                )
                local_fn = (
                    '4cker-output/{comparison}/{kind}_k{k}/'
                    '{bait}_{kind}_colorized_differential.bigbed'.format(**locals())
                )

                bait_bed_view.add_tracks(
                    Track(
                        name=sanitize(os.path.basename(local_fn)),
                        tracktype='bigBed 9',
                        itemRgb='on',
                        local_fn=local_fn,
                        short_label=os.path.basename(local_fn),
                        long_label=os.path.basename(local_fn),
                        subgroups=bait_subgroups,
                        color=color_for_track(local_fn),
                    )
                )
                continue

            agg_key = (comparison, kind, bait, treatment)
            if agg_key not in aggregates:
                aggregates[agg_key] = AggregateTrack(
                    name=''.join(agg_key),
                    tracktype='bigWig',
                    short_label=' '.join(agg_key),
                    long_label=' '.join(agg_key),
                    aggregate='transparentOverlay',
                    viewLimits='0:150000',
                    maxHeightPixels='8:50:127',
                    visibility='full',
                )
                supertrack.add_track(aggregates[agg_key])
            aggregate = aggregates[agg_key]

            local_fn = (
                '4cker-output/{comparison}/{kind}_k{k}/'
                '{bait}_{treatment}_{kind}_highinter.bigbed'.format(**locals())
            )

            bait_bed_view.add_tracks(
                Track(
                    name=sanitize(os.path.basename(local_fn)),
                    tracktype='bigBed 3',
                    local_fn=local_fn,
                    short_label=os.path.basename(local_fn),
                    long_label=os.path.basename(local_fn),
                    subgroups=bait_subgroups,
                    color=color_for_track(local_fn),
                )
            )

            for sample in vals[treatment]:
                for inter in ['noninter', 'lowinter', 'highinter', 'any']:

                    sample_subgroups = dict(
                        bait=sanitize(bait),
                        treatment=sanitize(treatment),
                        sample=sanitize(sample),
                        kind=sanitize(kind),
                        inter=sanitize(inter),
                    )

                    for key, val in sample_subgroups.items():
                        sample_subgroup_values[sanitize(key)].update(
                            [sanitize(val, False)])

                    if inter == 'any':
                        local_fn = (
                            '4cker-output/{comparison}/{kind}_k{k}/'
                            '{sample}_{kind}_norm_counts.bigwig'.format(**locals())
                        )


                        aggregate.add_subtrack(
                            Track(
                                name=os.path.basename(local_fn).replace('_', '').replace('.', ''),
                                tracktype='bigWig',
                                local_fn=local_fn,
                                short_label=os.path.basename(local_fn),
                                long_label=os.path.basename(local_fn),
                                subgroups=sample_subgroups,
                                color=color_for_track(local_fn),
                            )
                        )
                        sample_signal_view.add_tracks(
                            Track(
                                name=os.path.basename(local_fn).replace('_', '').replace('.', '') + 'smp',
                                tracktype='bigWig',
                                local_fn=local_fn,
                                short_label=os.path.basename(local_fn),
                                long_label=os.path.basename(local_fn),
                                subgroups=sample_subgroups,
                                color=color_for_track(local_fn),
                            )
                        )

                        continue
                    local_fn = (
                        '4cker-output/{comparison}/{kind}_k{k}/'
                        '{sample}_{kind}_{inter}.bigbed'.format(**locals())
                    )

                    sample_bed_view.add_tracks(
                        Track(
                            name=os.path.basename(local_fn).replace('_', '').replace('.', ''),
                            tracktype='bigBed 3',
                            local_fn=local_fn,
                            visibility='dense',
                            short_label=os.path.basename(local_fn),
                            long_label=os.path.basename(local_fn),
                            subgroups=sample_subgroups,
                            color=color_for_track(local_fn),
                        )
                    )



sample_subgroups = []
for k, v in sample_subgroup_values.items():
    sample_subgroups.append(
        SubGroupDefinition(
            name=sanitize(k),
            label=sanitize(k, strict=False),
            mapping={sanitize(i): sanitize(i, False) for i in list(v)}
        )
    )
sample_composite.add_subgroups(sample_subgroups)
sample_composite.add_params(
    dimensions='dimX=bait dimY=kind dimA=sample dimB=inter',
    sortOrder='bait=+ sample=+',
    filterComposite='dimA dimB',
)


bait_subgroups = []
for k, v in bait_subgroup_values.items():
    bait_subgroups.append(
        SubGroupDefinition(
            name=sanitize(k),
            label=sanitize(k, strict=False),
            mapping={sanitize(i): sanitize(i, False) for i in list(v)}
        )
    )
bait_composite.add_subgroups(bait_subgroups)
bait_composite.add_params(
    dimensions='dimX=bait dimY=kind dimA=treatment',
    sortOrder='bait=+',
    filterComposite='dimA',
)



f = hub.render()
linknames = upload.upload_hub(
    host=config['trackhub']['upload']['host'],
    user=config['trackhub']['upload']['user'],
    hub=hub
)
