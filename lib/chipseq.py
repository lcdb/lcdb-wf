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
