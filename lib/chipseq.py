"""
Helpers for ChIP-seq.
"""

# Example config for reference
# __example_config__ = {
#     'peak_calling': {
#         [
#             {
#                 'label': 'rep1',
#                 'algorithm': 'macs2',
#                 'input': ['input_1'],
#                 'ip': ['ip_1'],
#                 'extra': '--gs dm',
#             },
#             {
#                 'label': 'rep2',
#                 'algorithm': 'macs2',
#                 'input': ['input_2'],
#                 'ip': ['ip_2'],
#                 'extra': '--gs dm',
#             },
#
#         ]
#     }
# }


def peak_calling_dict(config, algorithm=None):
    """
    Returns a dictionary of peak-calling runs from the config.

    Parameters
    ----------
    config : dict

    algorithm : None
        If `algorithm` is None, dictionary is keyed by (label, algorithm).
        Otherwise, only the runs for `algorithm` are returned, keyed by
        `label`.
    """
    d = {}

    if 'chipseq' not in config:
        return d

    if config['chipseq'] is None:
        return d

    peaks_blocks = config['chipseq'].get('peak_calling', [])
    if not peaks_blocks:
        return d

    for block in peaks_blocks:
        key = (block['label'], block['algorithm'])
        if algorithm:
            if key[1] != algorithm:
                continue
            key = key[0]
        if key in d:
            raise ValueError("peak calling run '{0}' already defined".format(key))

        # If metadata key has been provided, then use that to populate the
        # block as default values.
        metadata = config['references'][config['organism']][config['aligner']['tag']].get('metadata', {})
        block.update(metadata)

        d[key] = block
    return d


def block_for_run(config, label, algorithm):
    """
    Returns the block for the (label, algorithm) run.

    Parameters
    ----------
    config : dict

    label, algorithm : str
    """
    return peak_calling_dict(config)[(label, algorithm)]


def samples_for_run(config, label, algorithm, treatment):
    """
    Returns the sample names configured for a particular peak-calling run

    Parameters
    ----------
    config : dict

    label, algorithm :
        Used as keys into peak_calling_dict()

    treatment : ip | input
    """
    return block_for_run(config, label, algorithm)[treatment]


def merged_input_for_ip(sampletable, merged_ip):
    """
    Returns the merged input label for a merged IP label.

    This is primarily used for the `fingerprint` rule, where we collect all the
    available input BAMs together.

    Parameters
    ----------

    sampletable : pandas.DataFrame

    merged_ip : str
        Label of IP to use, must be present in the `label` column of the
        sampletable.

    Examples
    --------
    This should make more sense if we have an example to work with.....

    Samples ip1 and ip2 are technical replicates. They are from a different
    experiment than ip3 and input3, hence their different
    biological_material.

    The way we know that input1 should be paired with ip1 and ip2 is because
    it shares the same biological material.

    Compare input1 and input9. They are not technical replicates (since they
    do not share the same `label`) but they are biological replicates because
    they share the same biological material.

    >>> from io import StringIO
    >>> import pandas as pd
    >>> df = pd.read_csv(StringIO('''
    ... samplename  antibody   biological_material  label
    ... ip1         gaf        s2cell-1             s2cell-gaf-1
    ... ip2         gaf        s2cell-1             s2cell-gaf-1
    ... ip3         ctcf       s2cell-2             s2cell-ctcf-1
    ... input1      input      s2cell-1             s2cell-input-1
    ... input3      input      s2cell-2             s2cell-input-3
    ... input9      input      s2cell-1             s2cell-input-1'''),
    ... sep='\s+')


    >>> merged_input_for_ip(df, 's2cell-gaf-1')
    ['s2cell-input-1']

    >>> merged_input_for_ip(df, 's2cell-ctcf-1')
    ['s2cell-input-3']

    """
    biomaterial = sampletable.loc[sampletable['label'] == merged_ip, 'biological_material']

    # This is a double-check that the merged IP all comes from the same
    # biological material.
    #
    # Note that we're printing here in addition to raising a ValueError because
    # when this is used in an input function to a snakemake rule, the traceback
    # is hidden
    #
    if biomaterial.nunique() != 1:
        print('Expected a single biomaterial, found: {0}'.format(biomaterial))
        raise ValueError

    # Now we find all the inputs from that biological material -- this is
    # defined as everything sharing the biological material that has "input" as
    # its antibody. They should also
    biomaterial = biomaterial.values[0]
    input_label = sampletable.loc[
        (sampletable['biological_material'] == biomaterial) &
        (sampletable['antibody'] == 'input'),
        'label'
    ]

    return sorted(input_label.unique())


def detect_peak_format(fn):
    """
    Figure out if a BED file is narrowPeak or broadPeak.

    Returns None if undetermined.

    This is useful for figuring out which autoSql file we should use or which
    bigBed 6+4 or bigBed 6+3 format to use.
    """
    line = open(fn).readline().strip()
    toks = line.split('\t')
    if len(toks) == 10:
        return 'narrowPeak'
    if len(toks) == 9:
        return 'broadPeak'
