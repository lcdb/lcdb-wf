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
from . import utils

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

        self.config = common.load_config(
            common.resolve_config(config, workdir))

        # Read the config file and extract all sort of useful bits. This mostly
        # uses the `common` module to handle the details.
        self.config['references_dir'] = common.get_references_dir(self.config)
        self.samples, self.sampletable = common.get_sampletable(self.config)
        self.refdict, self.conversion_kwargs = common.references_dict(self.config)
        self.organism = self.config['organism']

        # Patterns have one of three structures:
        #
        # Structure 1 (simple):
        #
        #    salmon:
        #      pattern: 'data/rnaseq_samples/{sample}/{sample}.salmon/{sample}_quant.sf'
        #      description: 'Transcripts quantification using Salmon'
        #
        # Structure 2 (nested 1 deep):
        #
        #    rseqc:
        #      bam_stat:
        #        pattern: 'data/rnaseq_samples/{sample}/rseqc/{sample}_bam_stat.txt'
        #        description: 'RNAseq quality control analysis'
        #      infer_experiment:
        #        pattern: 'data/rnaseq_samples/{sample}/rseqc/{sample}_infer_experiment.txt'
        #        description: 'Infer layout and strandedness of experiment'
        #
        #
        # Structure 3 (nested 2 deep):
        #
        #    patterns_by_peaks:
        #      peaks:
        #        macs2:
        #          pattern: 'data/chipseq_peaks/macs2/{macs2_run}/peaks.bed'
        #          description: 'BED file of peaks from macs2 peak-caller'
        #
        #
        # Our job here is to create three objects:
        #
        # patterns = {
        #   'salmon': 'data/rnaseq_samples/{sample}/{sample}.salmon/{sample}_quant.sf',
        #   'rseqc': {
        #       'bam_stat': 'data/rnaseq_samples/{sample}/rseqc/{sample}_bam_stat.txt',
        #       'infer_experiment': 'data/rnaseq_samples/{sample}/rseqc/{sample}_infer_experiment.txt',
        #   },
        #   'patterns_by_peaks': {
        #     'peaks': {
        #       'macs2': 'data/chipseq_peaks/macs2/{macs2_run}/peaks.bed'
        #     }
        #   }
        # }
        #
        # rst_files = {
        #   'salmon': 'reports/data/rnaseq_samples/sample/sample.salmon/sample_quant.sf.rst',
        #   'rseqc': {
        #       'bam_stat': 'reports/data/rnaseq_samples/sample/rseqc/sample_bam_stat.txt',
        #       'infer_experiment': 'reports/data/rnaseq_samples/sample/rseqc/sample_infer_experiment.txt.rst',
        #   },
        #   'patterns_by_peaks': {
        #     'peaks': {
        #       'macs2': 'reports/data/chipseq_peaks/macs2/macs2_run/peaks.bed.rst'
        #     }
        # }
        #
        # descriptions = {
        #   'salmon': 'Transcripts quantification using Salmon'
        #   'rseqc': {
        #       'bam_stat': 'RNAseq quality control analysis',
        #       'infer_experiment': 'Infer layout and strandedness of experiment',
        #   },
        #   'patterns_by_peaks': {
        #     'peaks': {
        #       'macs2': 'BED file of peaks from macs2 peak-caller'
        #     }
        #   }
        # }
        #
        #


        loaded_patterns = yaml.load(open(patterns), Loader=yaml.FullLoader)

        self._loaded_patterns = loaded_patterns
        self.patterns = utils.extract_nested(loaded_patterns, "pattern")
        self.descriptions = utils.extract_nested(loaded_patterns, "description")
        self.rst_files = utils.map_nested_dicts(self.patterns, utils.pattern_to_rst_file)

        self.is_paired = helpers.detect_layout(self.sampletable) == "PE"
        if self.is_paired:
            self.n = [1, 2]
        else:
            self.n = [1]


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

        # The merged bigwigs have different placeholders and therefore must be
        # filled separately. They also only should be included if
        # merged_bigwigs have been configured.
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

            self.rst_files.update(self.rst_files.pop('patterns_by_aggregate'))
            self.descriptions.update(self.descriptions.pop('patterns_by_aggregate'))


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
        self.descriptions.update(self.descriptions['patterns_by_sample'])
        self.rst_files.update(self.rst_files['patterns_by_sample'])

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

            self.rst_files.update(self.rst_files.pop('patterns_by_aggregate'))
            self.descriptions.update(self.descriptions.pop('patterns_by_aggregate'))

        # Then the peaks...
        #
        # Note: when adding support for new peak callers, add them here.
        PEAK_CALLERS = ['macs2', 'spp', 'sicer']

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
        self.descriptions.update(self.descriptions['patterns_by_peaks'])
        self.rst_files.update(self.rst_files['patterns_by_peaks'])
