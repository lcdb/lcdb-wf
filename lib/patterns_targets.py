import os
import collections
import yaml
from . import common
from . import chipseq
from lcdblib.snakemake import helpers, aligners
from lcdblib.utils import utils

HERE = os.path.abspath(os.path.dirname(__file__))


def update_recursive(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = update_recursive(d.get(k, {}), v)
        else:
            d[k] = v
    return d


class RNASeqConfig(object):
    def __init__(self, config):

        self.path = None
        if isinstance(config, str):
            self.path = config

        configdict, pth = common.resolve_config(config)

        # Called for its side-effect of checking that a references dir is specified
        common.get_references_dir(config)

        self.samples, self.sampletable = common.get_sampletable(config)
        self.refdict, self.conversion_kwargs = common.references_dict(config)

        self.assembly = configdict['assembly']

        self.sample_dir = configdict.get('sample_dir', 'samples')
        self.agg_dir = configdict.get('aggregation_dir', 'aggregation')

        self.patterns = yaml.load(open(os.path.join(HERE, 'rnaseq_patterns.yaml')))
        self.fill = dict(sample=self.samples, sample_dir=self.sample_dir, agg_dir=self.agg_dir)
        self.targets = helpers.fill_patterns(self.patterns, self.fill)


class ChIPSeqConfig(object):
    def __init__(self, config):

        self.path = None
        if isinstance(config, str):
            self.path = config

        configdict, pth = common.resolve_config(config)

        # Called for its side-effect of checking that a references dir is specified
        common.get_references_dir(config)

        self.samples, self.sampletable = common.get_sampletable(config)
        self.refdict, self.conversion_kwargs = common.references_dict(config)

        self.assembly = configdict['assembly']

        self.sample_dir = configdict.get('sample_dir', 'samples')
        self.agg_dir = configdict.get('aggregation_dir', 'aggregation')
        self.merged_dir = configdict.get('merged_dir', 'merged')
        self.peak_calling = configdict.get('peaks_dir', 'chipseq')

        _patterns = yaml.load(open(os.path.join(HERE, 'chipseq_patterns.yaml')))

        self.patterns_by_sample = _patterns['patterns_by_sample']
        self.fill_by_sample = dict(sample=self.samples, sample_dir=self.sample_dir, agg_dir=self.agg_dir)
        self.fill_by_sample = dict(
            sample=self.samples.values, sample_dir=self.sample_dir, agg_dir=self.agg_dir,
            merged_dir=self.merged_dir, peak_calling=self.peak_calling,
            label=self.sampletable.label.values,
            ip_label=self.sampletable.label[self.sampletable.antibody != 'input'].values)
        self.targets_by_sample = helpers.fill_patterns(self.patterns_by_sample, self.fill_by_sample)

        self.patterns_by_peaks = _patterns['patterns_by_peaks']

        self.fill_by_peaks = dict(
            peak_calling=self.peak_calling,
            macs2_run=list(chipseq.peak_calling_dict(dict(configdict), algorithm='macs2').keys()),
            spp_run=list(chipseq.peak_calling_dict(dict(configdict), algorithm='spp').keys()),
            combination='zip',
        )

        self.targets_for_peaks = {}
        for pc in ['macs2', 'spp']:
            _peak_patterns = {}
            for k, v in self.patterns_by_peaks.items():
                _peak_patterns[k] = {pc: self.patterns_by_peaks[k][pc]}
            _fill = {
                'peak_calling': self.peak_calling,
                pc + '_run': list(chipseq.peak_calling_dict(dict(configdict), algorithm=pc).keys())}
            update_recursive(self.targets_for_peaks, helpers.fill_patterns(_peak_patterns, _fill))

        self.targets = {}
        self.targets.update(self.targets_by_sample)
        self.targets.update(self.targets_for_peaks)
        self.patterns = {}
        self.patterns.update(self.patterns_by_sample)
        self.patterns.update(self.patterns_by_peaks)

