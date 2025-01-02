import yaml
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import math


def lcdbwf_samplename(x):
    """
    Processes sample names by removing specific substrings and extracting relevant parts.

    Parameters
    ----------
    x : list of str
        A list of sample name strings.

    Returns
    -------
    list of str
        A list of processed sample names.
    """

    processed_names = []
    for name in x:
        # Remove unnecessary substrings and split to obtain the core sample name
        name = name.replace('data/rnaseq_samples/', '').replace('.cutadapt.bam', '')
        name_parts = name.split('/')
        processed_names.append(name_parts[0])
    return processed_names


def process_featurecounts(fc_path, sampletable_path, output_path, sample_func=lcdbwf_samplename, subset_counts=False):
    """
    Processes a feature counts table and matches it with a sample table.

    Parameters
    ----------
    fc_path : str
        Path to the feature counts file.
    sampletable_path : str
        Path to the sample table TSV file.
    output_path : str
        Path to save the processed counts DataFrame.
    sample_func : function, optional
        Function to process sample names from counts data. Default is `lcdbwf_samplename`.
    subset_counts : bool, optional
        If True, subsets counts to samples present in the sample table. Default is False.

    Returns
    -------
    None
    """

    # Load counts, set 'Geneid' as index, remove metadata columns, and process sample names
    m = pd.read_csv(fc_path, sep='\t', comment='#').set_index('Geneid').iloc[:, 5:]
    x = sample_func(m.columns.tolist())
    sampletable = pd.read_csv(sampletable_path, sep='\t')
    samplenames = sampletable.iloc[:, 0].tolist()

    # Identify samples mismatched between counts and sampletable
    counts_not_sampletable = set(x) - set(samplenames)
    sampletable_not_counts = set(samplenames) - set(x)

    # Handle mismatches and subset counts if specified
    if counts_not_sampletable:
        if not subset_counts:
            raise ValueError(
                f"The following samples are in the counts data but not in the sample table. "
                f"Consider using subset_counts=True to remove them from the counts: "
                f"{', '.join(counts_not_sampletable)}"
            )
        else:
            indices_to_keep = [i for i, name in enumerate(x) if name in samplenames]
            x = [name for name in x if name in samplenames]
            m = m.iloc[:, indices_to_keep]

    if sampletable_not_counts:
        raise ValueError(
            f"The following samples are in the sample table but not in the counts data. "
            f"Check sample_func? {', '.join(sampletable_not_counts)}"
        )

    # Align columns in m to sampletable order and save to output path
    m.columns = x
    m = m[samplenames]
    m.to_csv(output_path, sep='\t')


def geometric_mean(arr):
    """
    Calculates the geometric mean of non-zero elements of an array.

    Parameters
    ----------
    arr : Array of values.

    Returns
    -------
    float
        Geometric mean of non-zero elements in the array.
    """
    return np.exp(np.mean(np.log(arr)))


def estimate_sizefactors(counts_path, sizefactors_path, sigfigs = 6):
    """
    Estimates size factors for normalization using the median of ratios of counts to geometric means.

    Parameters
    ----------
    counts_path : str
        Path to the processed counts TSV file.
    sizefactors_path : str
        Path to save the estimated size factors as a TSV file.

    Returns
    -------
    None
    """

    # Load counts, remove rows with zeros, and compute geometric means for each gene
    counts_df = pd.read_csv(counts_path, sep='\t', index_col=0)
    counts_df = counts_df[counts_df > 0].dropna()
    geometric_means = counts_df.apply(geometric_mean, axis=1)

    # Compute size factors by taking median of ratios of counts to geometric means
    ratios = counts_df.div(geometric_means, axis=0)
    sizefactors_df = ratios.apply(np.median, axis=0)

    # Format size factors as DataFrame, reorder columns, and save
    sizefactors_df = pd.DataFrame(sizefactors_df, columns=['sizeFactors'])
    sizefactors_df['samplename'] = sizefactors_df.index
    sizefactors_df = sizefactors_df[['samplename', 'sizeFactors']].reset_index(drop=True)
    # Round size factors to 6 decimal places
    sizefactors_df['sizeFactors'] = sizefactors_df['sizeFactors'].apply(lambda x: round(x, 6))
    sizefactors_df.to_csv(sizefactors_path, sep='\t', header=True, index=False)


def normalize_user_region_counts_list(sizefactors_path, counts):
    """
    Normalize counts for specified regions using size factors.

    Parameters
    ----------
    sizefactors_path : str
        Path to the size factors TSV file.
    counts : list of str
        List of paths to counts files for different regions.

    Returns
    -------
    defaultdict
        Nested dictionary with normalized counts by sample and region.
    """

    sizefactors_df = pd.read_csv(sizefactors_path, sep='\t', index_col='samplename')
    normalized_counts = defaultdict(dict)

    for file_path in counts:
        parts = file_path.strip().split(os.sep)
        sample = parts[-2]
        region = parts[-1].replace('.counts.txt', '')

        with open(file_path) as infile:
            count = int(infile.read())

        size_factor = sizefactors_df.loc[sample, 'sizeFactors']
        # Normalize and round counts to 3 decimal places
        normalized_counts[sample][region] = round(float(count) / float(size_factor), 3)

    return normalized_counts


def make_regions_counts_yaml(sizefactors_path, counts_list, regions, yaml_path):
    """
    Generates a YAML file with normalized counts for specified regions for visualization in MultiQC.

    Parameters
    ----------
    sizefactors_path : str
        Path to the size factors TSV file.
    counts : list of str
        List of paths to counts files for different regions.
    regions : list of str
        Configured region names
    yaml_path : str
        Path to save the generated YAML file.

    Returns
    -------
    None
    """
    normalized_counts = normalize_user_region_counts_list(sizefactors_path, counts_list)
    data = {sample: regions_data for sample, regions_data in normalized_counts.items()}
    y = {
        'id': 'expression_barchart',
        'section_name': 'Normalized Counts in User Defined Regions',
        'description': f"The user-defined regions are: {', '.join(regions)}",
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'user_regions_counts_plot',
            'title': f"Normalized counts in {', '.join(regions)}",
            'ylab': 'Normalized Counts',
            'xlab': 'Samples'
        },
        'data': data
    }

    with open(yaml_path, 'w') as outfile:
        yaml.dump(y, outfile, default_flow_style=False)


def normalize_featurecounts(counts_path, sizefactors_path, output_path):
    """
    Normalizes the processed featurecounts table by size factors for each sample.

    Parameters
    ----------
    counts_path : str
        Path to the counts TSV file.
    sizefactors_path : str
        Path to the size factors TSV file.
    output_path : str
        Path to save the normalized counts TSV file.

    Returns
    -------
    None
    """

    # Read counts and size factors
    counts_df = pd.read_csv(counts_path, sep='\t', index_col=0)
    sizefactors_df = pd.read_csv(sizefactors_path, sep='\t', index_col='samplename')

    # Normalize each sample's counts by its size factor
    for sample in counts_df.columns:
        if sample in sizefactors_df.index:
            counts_df[sample] = round(counts_df[sample] / sizefactors_df.loc[sample, 'sizeFactors'], 3)
        else:
            raise ValueError(f"Size factor for sample '{sample}' not found in sizefactors file: '{sizefactors_path}'.")

    # Save the normalized counts table
    counts_df.to_csv(output_path, sep='\t')

def get_regions(config, c, final_targets):
    """
    Reads regions and their start and stop positions from config.
    Saves the user configured loci to a dictionary and uses the dictionary
    to make a list of region names. Updates final targets to include the expression
    multiqc module.

    Parameters
    __________
    config : dictionary
        config object from config/config.yaml
    c : what is c
        object of dictionary of target paths
    final_targets : list of strings
        list of paths to targets

    Returns
    _______
    Tuple of
    REGIONS : list of str
        region names defined in config
    REGION_DICT : dict of str
        region names and chrom:start-stop
    final_targets : list of str
        target output files
    """

    final_targets.append(c.targets['expression_mqc_yaml'])
    REGION_DICTS = []
    for name, coord in config['expression_barchart_multiqc'].items():
        chrom, positions = coord.strip().split(':')
        start, end = positions.strip().split('-')
        start = int(start)
        end = int(end)
        REGION_DICTS.append({'name': name, 'chrom': chrom, 'start': start, 'end': end})

    REGIONS = [region['name'] for region in REGION_DICTS]

    return (REGIONS, REGION_DICTS, final_targets)

