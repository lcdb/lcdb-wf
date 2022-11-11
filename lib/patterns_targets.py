"""
This module handles the reading and filling-in of patterns. It can be used from
within Snakefiles or in downstream (figure-making) scripts.
"""

import os
import collections
import yaml
from . import common
from . import chipseq
from . import helpers

HERE = os.path.abspath(os.path.dirname(__file__))

# Note: when adding support for new peak callers, add them here.
PEAK_CALLERS = ['macs2', 'spp', 'sicer']


def update_recursive(d, u):
    """
    Update dictionary `d` with items in dictionary `u`, recursively
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_recursive(d.get(k, {}), v)
        else:
            d[k] = v
    return d


class SeqConfig(object):
    def __init__(self, config, patterns, workdir=None):
        """
        This class takes care of common tasks related to config and patterns
        files (reading the sampletable, etc) but is intended to be subclassed.

        Parameters
        ----------
        config : str or dict

        patterns : str
            Path to patterns YAML file

        workdir : str
            Config, patterns, and all paths in `config` should be interpreted
            as relative to `workdir`
        """
        self.path = None
        self.workdir = '.'
        if workdir is not None:
            config = os.path.join(workdir, config)
            patterns = os.path.join(workdir, patterns)
            self.workdir = workdir

        if isinstance(config, str):
            self.path = config

        self.config = common.load_config(
            common.resolve_config(config, workdir))

        stranded = self.config.get('stranded', None)
        self.stranded = None
        if stranded:
            if stranded in ('unstranded'):
                self.stranded = 'unstranded'
            elif stranded in ('fr-firststrand', 'ISR', 'SR', 'reverse'):
                self.stranded = 'fr-firststrand'
            elif stranded in ('fr-secondstrand', 'ISF', 'SF', 'forward'):
                self.stranded = 'fr-secondstrand'

        # Read the config file and extract all sort of useful bits. This mostly
        # uses the `common` module to handle the details.
        self.config['references_dir'] = common.get_references_dir(self.config)
        self.samples, self.sampletable = common.get_sampletable(self.config)
        self.refdict, self.conversion_kwargs = common.references_dict(self.config)
        self.organism = self.config['organism']
        self.patterns = yaml.load(open(patterns), Loader=yaml.FullLoader)
        self.is_paired = helpers.detect_layout(self.sampletable) == 'PE'
        if self.is_paired:
            self.n = [1, 2]
        else:
            self.n = [1]

        helpers.preflight(self.config)

class RNASeqConfig(SeqConfig):
    def __init__(self, config, patterns, workdir=None):
        """
        Config object specific to RNA-seq workflows.

        Fills in patterns to create targets by handling the by-sample and
        by-aggregate sections separately.

        Parameters
        ----------

        config : dict

        patterns : str
            Path to patterns YAML file

        workdir : str
            Config, patterns, and all paths in `config` should be interpreted
            as relative to `workdir`
        """
        SeqConfig.__init__(self, config, patterns, workdir)

        self.fill = dict(sample=self.samples, n=self.n)
        self.patterns_by_aggregation = self.patterns.pop('patterns_by_aggregate', None)
        self.targets = helpers.fill_patterns(self.patterns, self.fill, zip)

        # Then the aggregation
        if self.patterns_by_aggregation is not None and 'merged_bigwigs' in self.config:
            self.fill_by_aggregation = dict(
                merged_bigwig_label=self.config['merged_bigwigs'].keys(),
            )
            self.targets_by_aggregation = helpers.fill_patterns(
                self.patterns_by_aggregation,
                self.fill_by_aggregation
            )
            self.targets.update(self.targets_by_aggregation)
            self.patterns.update(self.patterns_by_aggregation)

        helpers.rnaseq_preflight(self)


class ChIPSeqConfig(SeqConfig):
    def __init__(self, config, patterns, workdir=None):
        """
        Config object specific to ChIP-seq workflows.

        Fills in patterns to create targets by handling the by-sample, by-peak,
        and by-aggregate sections separately.

        Parameters
        ----------

        config : dict

        patterns : str
            Path to patterns YAML file

        workdir : str
            Config, patterns, and all paths in `config` should be interpreted
            as relative to `workdir`
        """
        SeqConfig.__init__(self, config, patterns, workdir)

        self.targets = {}

        # For ChIP-seq, the structure of the patterns is quite different for
        # samples than it is for peaks. For example, the peaks do not have any
        # sample info in the filenames but aggregate possibly many different samples
        #
        # So construct them separately, and then later update self.patterns and
        # self.targets.
        #
        # The averaged bigwigs are also aggregated, but in a different way.
        # They will be handled separately.
        #
        # First, the samples...
        self.patterns_by_sample = self.patterns['patterns_by_sample']
        self.fill_by_sample = dict(
            n=self.n,
            sample=self.samples.values,
            label=self.sampletable.label.values,
            ip_label=self.sampletable.label[
                self.sampletable.antibody != 'input'].values
        )
        self.targets_by_sample = helpers.fill_patterns(
            self.patterns_by_sample, self.fill_by_sample)

        self.targets.update(self.targets_by_sample)
        self.patterns.update(self.patterns_by_sample)

        # Then the aggregation...
        self.patterns_by_aggregation = self.patterns.pop('patterns_by_aggregate', None)
        if self.patterns_by_aggregation is not None and 'merged_bigwigs' in self.config:
            self.fill_by_aggregation = dict(
                merged_bigwig_label=self.config['merged_bigwigs'].keys(),
            )
            self.targets_by_aggregation = helpers.fill_patterns(
                self.patterns_by_aggregation,
                self.fill_by_aggregation
            )
            self.targets.update(self.targets_by_aggregation)
            self.patterns.update(self.patterns_by_aggregation)

        # Then the peaks...
        #

        self.patterns_by_peaks = self.patterns['patterns_by_peaks']
        self.targets_for_peaks = {}

        # We need to fill in just those peak-calling runs that are specified
        # for each peak-caller. For reference, here's an example
        # `patterns_by_peaks` from the YAML:
        #
        #        peaks:
        #           macs2: '{peak_calling}/macs2/{macs2_run}/peaks.bed'
        #           spp: '{peak_calling}/spp/{spp_run}/peaks.bed'
        #        bigbed:
        #            macs2: '{peak_calling}/macs2/{macs2_run}/peaks.bigbed'
        #            spp: '{peak_calling}/spp/{spp_run}/peaks.bigbed'


        # Also note that the snakefile's all rule uses
        # utils.flatten(c.targets['peaks']), but in the case where no
        # peak-calling runs are specified these should be initialized,
        # otherwise we'll get a KeyError.
        self.targets['peaks'] = []
        self.targets['bigbed'] = []

        for pc in PEAK_CALLERS:
            # Extract out just the subset of `patterns_by_peaks` for this
            # peak-caller e.g., from the example above, if pc='macs2' this
            # would only be:
            #
            #   peaks:
            #      macs2: '{peak_calling}/macs2/{macs2_run}/peaks.bed'
            #   bigbed:
            #       macs2: '{peak_calling}/macs2/{macs2_run}/peaks.bigbed'
            #
            _peak_patterns = {
                k: {pc: v[pc]} for k, v in self.patterns_by_peaks.items()
            }


            # Fix for issue #166, which was caused by commit 8a211122:
            #
            # If no runs for the peak-caller are configured, this will be
            # empty and we should continue on.
            peaks_to_fill = list(chipseq.peak_calling_dict(self.config, algorithm=pc).keys())

            if not peaks_to_fill:
                continue

            _fill = {pc + '_run': peaks_to_fill}

            # The trick here is the recursive updating of targets_for_peaks.
            # We're adding the filled-in runs of each peak caller to the
            # targets as they're built.
            update_recursive(
                self.targets_for_peaks,
                helpers.fill_patterns(_peak_patterns, _fill)
            )


        self.targets.update(self.targets_for_peaks)
        self.patterns.update(self.patterns_by_peaks)

        helpers.chipseq_preflight(self)
