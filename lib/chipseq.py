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
    for block in config['chipseq']['peak_calling']:
        key = (block['label'], block['algorithm'])
        if algorithm:
            if key[1] != algorithm:
                continue
            key = key[0]
        if key in d:
            raise ValueError("peak calling run '{0}' already defined".format(key))
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

    Parameters
    ----------

    sampletable : pandas.DataFrame

    merged_ip : str
        Label of IP to use, must be present in the `label` column of the
        sampletable.
    """
    biomaterial = sampletable.loc[sampletable['label'] == merged_ip, 'biological_material']

    # note that we're printing here in addition to raising a ValueError because
    # when this is used in an input function to a snakemake rule, the traceback
    # is hidden
    if biomaterial.nunique() != 1:
        print('Expected a single biomaterial, found: {0}'.format(biomaterial))
        raise ValueError
    biomaterial = biomaterial.values[0]
    input_label = sampletable.loc[
        (sampletable['biological_material'] == biomaterial) &
        (sampletable['antibody'] == 'input'),
        'label'
    ]

    if input_label.nunique() != 1:
        print("Expected a single input label, found: {0}".format(input_label))

    return input_label.values[0]


def inputs_for_ip(sampletable, ip):
    """
    Returns all inputs with the same biological material as the provided IP
    sample.

    This can be useful for creating the "chipseq:peak_calling" section of the
    config.yaml.

    Differs from `merged_input_for_ip` in that here we return *all* inputs on
    a sample level (first column of the sampletable) while that function
    expects a 1:1 mapping between *merged* samples (the "label" column of the
    sampletable).
    """
    biomaterial = sampletable.loc[sampletable[sampletable.columns[0]] == ip, 'biomaterial']
    if biomaterial.nunique() != 1:
        print('Expected a single biomaterial, found: {0}'.format(biomaterial))
        raise ValueError
    biomaterial = biomaterial.values[0]
    inputs = sampletable.loc[
        (sampletable['biological_material'] == biomaterial) &
        (sampletable['antibody'] == 'input'),
        sampletable.columns[0]
    ]
    return list(inputs)
