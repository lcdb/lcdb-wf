"""
This module handles the reading and filling-in of patterns. It can be used from
within Snakefiles or in downstream (figure-making) scripts.

A "patterns" YAML file organizes filename patterns, like
"{sample_dir}/{sample}/{sample}.unique.bam", used by Snakemake rules. There are
several advantages of this over writing the patterns directly in each rule:

- Patterns can be filled in automatically by a sample table (in fact, that's
  what's happening in the classes here), which means that the entire workflow
  is largely configured by a TSV

- The filled-in patterns (which we call "targets") can be used for the first
  `all` rule of a Snakefile and to easily disable large parts of the workflow
  by commenting them out in the `all` rule's input.

- Aggregation rules can use filled-in targets just by specifying the keys. E.g,
  `targets['bam']` instead of
  `expand('{sample_dir}/{sample}/{sample}.cutdapt.bam', sample=samples)`

- It is easier to customize locations and directory structure when all filename
  patterns are in the same place, rather than scrolling around a Snakefile.

- The organization allows us to more easily find what files are stored where.
  Rather than scroll through the entire workflow trying to find where the
  postitive-strand bigwigs are, we can either access the pattern with
  `patterns['bigwig']['pos']` or get the list of all plus-strand bigwigs across
  all samples with `targets['bigwig']['pos']`.

- Keeping the patterns outside the Snakefile allows us to re-use the patterns
  and their organization for downstream work (like figure-making).

"""

import os
import collections
import yaml
from . import common
from . import chipseq
from lcdblib.snakemake import helpers

HERE = os.path.abspath(os.path.dirname(__file__))


def update_recursive(d, u):
    """
    Update dictionary `d` with items in dictionary `u`, recursively
    """
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
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

        self.config = common.resolve_config(config, workdir)

        # Read the config file and extract all sort of useful bits. This mostly
        # uses the `common` module to handle the details.
        self.config['references_dir'] = common.get_references_dir(self.config)
        self.samples, self.sampletable = common.get_sampletable(self.config)
        self.refdict, self.conversion_kwargs = common.references_dict(self.config)
        self.assembly = self.config['assembly']
        self.patterns = yaml.load(open(patterns))


class RNASeqConfig(SeqConfig):
    def __init__(self, config, patterns, workdir=None):
        """
        Config object specific to RNA-seq workflows.
        """
        SeqConfig.__init__(self, config, patterns, workdir)
        self.sample_dir = self.config.get('sample_dir', 'samples')
        self.agg_dir = self.config.get('aggregation_dir', 'aggregation')
        self.fill = dict(sample=self.samples, sample_dir=self.sample_dir,
                         agg_dir=self.agg_dir)
        self.targets = helpers.fill_patterns(self.patterns, self.fill)


class ChIPSeqConfig(SeqConfig):
    def __init__(self, config, patterns, workdir=None):
        """
        Config object specific to ChIP-seq workflows.
        """
        SeqConfig.__init__(self, config, patterns, workdir)
        self.sample_dir = self.config.get('sample_dir', 'samples')
        self.agg_dir = self.config.get('aggregation_dir', 'aggregation')
        self.merged_dir = self.config.get('merged_dir', 'merged')
        self.peak_calling = self.config.get('peaks_dir', 'chipseq')

        # For ChIP-seq, the structure of the patterns is quite different for
        # samples vs for peaks. For example, the peaks do not have any sample
        # info in the filenames and are sort of at an aggregated level (since
        # multiple samples can make it into the same peak-calling run). So we
        # construct them separately, and then later update self.patterns and
        # self.targets.
        #
        # First, the samples
        self.patterns_by_sample = self.patterns['patterns_by_sample']
        self.fill_by_sample = dict(
            sample=self.samples.values,
            sample_dir=self.sample_dir,
            agg_dir=self.agg_dir,
            merged_dir=self.merged_dir,
            peak_calling=self.peak_calling,
            label=self.sampletable.label.values,
            ip_label=self.sampletable.label[
                self.sampletable.antibody != 'input'].values
        )
        self.targets_by_sample = helpers.fill_patterns(
            self.patterns_by_sample, self.fill_by_sample)

        # Then the peaks
        #
        # Note: when adding support for new peak callers, add them here.
        PEAK_CALLERS = ['macs2', 'macs2_extended', 'spp', 'sicer']

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

        for pc in PEAK_CALLERS:
            # extract out just the subset of `patterns_by_peaks` for this
            # peak-caller e.g.,
            #
            #   peaks:
            #      macs2: '{peak_calling}/macs2/{macs2_run}/peaks.bed'
            #   bigbed:
            #       macs2: '{peak_calling}/macs2/{macs2_run}/peaks.bigbed'
            #
            _peak_patterns = {
                k:{pc: v[pc]} for k, v in self.patterns_by_peaks.items()
            }


            _fill = {
                'peak_calling': self.peak_calling,
                pc + '_run': list(
                    chipseq.peak_calling_dict(self.config, algorithm=pc).keys())
            }

            # The trick here is the recursive updating of the
            # targets_for_peaks, so we're adding the filled-in runs of each
            # peak caller to the targets as they're built.
            update_recursive(
                self.targets_for_peaks,
                helpers.fill_patterns(_peak_patterns, _fill)
            )

        self.targets = {}
        self.targets.update(self.targets_by_sample)
        self.targets.update(self.targets_for_peaks)

        self.patterns = {}
        self.patterns.update(self.patterns_by_sample)
        self.patterns.update(self.patterns_by_peaks)
